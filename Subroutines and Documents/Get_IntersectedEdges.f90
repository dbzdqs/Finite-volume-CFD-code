!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: May., 15, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Get_IntersectedEdges(Dim,NC,Corn,Neib,X,Y,FrontEdges,States,Fronts,N_c,N_d,IntersectedEdges,IECount,TopEdgeRecoveryIsPossible) ! Implementation of Algorithm 2 in the paper
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,N_c,N_d,Corn,Neib,X,Y,FrontEdges,States,Fronts
Intent(Out)::IntersectedEdges,IECount,TopEdgeRecoveryIsPossible

Integer::Dim,NC,N_c,N_d,TEC,QEC,N_c_index,E,IECount,Ti,Ti1,I,J,Ni,corner1,corner2,Fronts,prvC,nxtC,oppC,index
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::IntersectedEdges
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::IntersectionOccur,isFrontEdge,TopEdgeRecoveryIsPossible,IsQuadDiagonal,isNotChevron,PointLaysOnLineSegment,IsDiagonalInsideQuad
!Real(8)::a,b,c
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
TopEdgeRecoveryIsPossible = .True.

if(N_c /= N_d) then
!Part 1:
	!__________________________ Finding Elements surrounding point N_c _________________________
	E = 0
	do I=1,NC
		if(Corn(I,1) /= Corn(I,2)) then
			do J=1,4
				if(Corn(I,J) == N_c) then
					E = I
					exit
				endif
			end do
		endif

		if(E /= 0) exit

	end do

	Call GetSurroundingElements(Dim,Corn,Neib,N_c,E,TElms,QElms,TEC,QEC)
!Part 2:
	!___________________________________ Finding triangle Ti ___________________________________
	do I=1,TEC

		if(.Not. TopEdgeRecoveryIsPossible) then
			
			print*,'>>>>> Unfortunately Top Edge Recovery is Impossible!!! <<<<<'
			exit

		endif

		do J=1,3
			if(Corn(TElms(I),J) == N_c) then 
				N_c_index = J
			endif
		end do
		!------------------------------- Finding other corners ----------------------------------
		if(N_c_index == 1) then
			corner1 = Corn(TElms(I),2)
			corner2 = Corn(TElms(I),3) 
		elseif(N_c_index == 2) then
			corner1 = Corn(TElms(I),1)
			corner2 = Corn(TElms(I),3)
		elseif(N_c_index == 3) then
			corner1 = Corn(TElms(I),1)
			corner2 = Corn(TElms(I),2)
		endif
!Part 3:
		if((corner1/=N_c .And. corner1/=N_d .And. PointLaysOnLineSegment(Dim,X,Y,N_c,N_d,corner1)) .Or. (corner2/=N_c .And. corner2/=N_d .And. PointLaysOnLineSegment(Dim,X,Y,N_c,N_d,corner2))) then
			
			print*,'Fail!!!'
			TopEdgeRecoveryIsPossible = .False.
			exit
		
		else

			!-------------------------- Situation in which intersection occurs ----------------------	
			if(IntersectionOccur(Dim,N_c,N_d,corner1,corner2,X,Y)) then
				Ti = TElms(I)
				!-------------------------- Finding Ei the edge opposite N_c ------------------------
				if(N_c_index == 1) then
					Ei(1) = Corn(Ti,2)
					Ei(2) = Corn(Ti,3)
				elseif(N_c_index == 2) then
					Ei(1) = Corn(Ti,1)
					Ei(2) = Corn(Ti,3)
				elseif(N_c_index == 3) then
					Ei(1) = Corn(Ti,1)
					Ei(2) = Corn(Ti,2)
				endif
				!-- Checking if Ei is not on the front in order to add it to list of intersecting edges--
				if(.Not. isFrontEdge(Dim,Fronts,FrontEdges,States,Ei)) then
					IECount = IECount + 1
					IntersectedEdges(IECount,1) = Ei(1)
					IntersectedEdges(IECount,2) = Ei(2)
		print *,'Tri:',Ti,'Intersected Edge:',IntersectedEdges(IECount,1),IntersectedEdges(IECount,2)
				else
					print*,'Fail!!!'
					TopEdgeRecoveryIsPossible = .False.
					exit
                endif
				!_____________________________ Adding other Edges intersecting Edge S ______________________
!Part 4:                
				do 	
					Call GetNextTriangle(Dim,NC,Corn,Ti,Ei,Ti1)
					
					if(Corn(Ti1,1)==N_d .Or. Corn(Ti1,2)==N_d .Or. Corn(Ti1,3)==N_d) then
						exit	
					endif

					Ti=Ti1

					do J=1,3
						if(Corn(Ti,J)/=Ei(1) .And. Corn(Ti,J)/=Ei(2)) then
							Ni = J
							print *,'Ni',Ni,':',Corn(Ti,J)
						endif
					end do
					
					if(Ni == 1) then
						if(IntersectionOccur(Dim,N_c,N_d,Corn(Ti,3),Corn(Ti,Ni),X,Y)) then
							Ei(1) = Corn(Ti,3)
							Ei(2) = Corn(Ti,Ni)
						else
							Ei(1) = Corn(Ti,Ni)
							Ei(2) = Corn(Ti,2)
						endif
					elseif(Ni == 2) then
						if(IntersectionOccur(Dim,N_c,N_d,Corn(Ti,1),Corn(Ti,Ni),X,Y)) then
							Ei(1) = Corn(Ti,1)
							Ei(2) = Corn(Ti,Ni)
						else
							Ei(1) = Corn(Ti,Ni)
							Ei(2) = Corn(Ti,3)
						endif
					elseif(Ni == 3) then
						if(IntersectionOccur(Dim,N_c,N_d,Corn(Ti,2),Corn(Ti,Ni),X,Y)) then
							Ei(1) = Corn(Ti,2)
							Ei(2) = Corn(Ti,Ni)
						else
							Ei(1) = Corn(Ti,1)
							Ei(2) = Corn(Ti,Ni)
						endif
					endif

					!-- Checking if Ei is not on the front in order to add it to list of intersecting edges--
					if(.Not. isFrontEdge(Dim,Fronts,FrontEdges,States,Ei)) then
						IECount = IECount + 1
						IntersectedEdges(IECount,1) = Ei(1)
						IntersectedEdges(IECount,2) = Ei(2)
		print *,'Tri:',Ti,'Intersected Edge:',IntersectedEdges(IECount,1),IntersectedEdges(IECount,2)
					else
						print*,'Fail!!!'
						TopEdgeRecoveryIsPossible = .False.
						exit
					endif
				end do
			endif

		endif

	end do
!Part 5:
	if(TopEdgeRecoveryIsPossible) then
		
		do I=1,QEC

			E = QElms(I)

			if(IsQuadDiagonal(Dim,Corn,N_c,N_d,E)) then

				if(isNotChevron(Dim,Corn,X,Y,E)) then

					print*,'Fail!!!'
					TopEdgeRecoveryIsPossible = .False.
					exit
				else

					if(IsDiagonalInsideQuad(Dim,Corn,X,Y,E,N_c,N_d)) then
						
						print*,'Fail!!!'
						TopEdgeRecoveryIsPossible = .False.
						exit

					endif

				endif

			else

				do J=1,4
					if(Corn(E,J) == N_c) then
						 index = J
						 exit
					endif
				end do

				if(index == 1) then
					prvC = Corn(E,4)
					nxtC = Corn(E,2)
					oppC = Corn(E,3)
				elseif(index == 2) then
					prvC = Corn(E,1)
					nxtC = Corn(E,3)
					oppC = Corn(E,4)
				elseif(index == 3) then
					prvC = Corn(E,2)
					nxtC = Corn(E,4)
					oppC = Corn(E,1)
				elseif(index == 4) then
					prvC = Corn(E,3)
			BAADtC = Corn(E,1)
					oppC = Corn(E,2)
				endif

				if(IntersectionOccur(Dim,N_c,N_d,prvC,oppC,X,Y)) then
					print*,'Fail!!!'
					TopEdgeRecoveryIsPossible = .False.
					exit
				endif

				if(IntersectionOccur(Dim,N_c,N_d,nxtC,oppC,X,Y)) then
					print*,'Fail!!!'
					TopEdgeRecoveryIsPossible = .False.
					exit
				endif

			endif

		end do

	endif

else
	
	TopEdgeRecoveryIsPossible = .False.

endif
!===========================================================================================
End Subroutine Get_IntersectedEdges
!*********************************************************************************************
