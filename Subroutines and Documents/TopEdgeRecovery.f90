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
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine TopEdgeRecovery(Dim,NC,N_c,N_d,Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible) ! Tested
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,N_c,N_d,X,Y,States,Fronts
Intent(InOut)::Corn,Neib,FrontEdges,TopEdgeRecoveryIsPossible

Integer::Dim,NC,Fronts,N_c,N_d,P1,P2,P3,P4,I,J,K,L,M,I1,I2,Iface1,Iface2,tri1,tri2,IECount
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::IntersectedEdges
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::IntersectionOccur,SwapIsSafe,TopEdgeRecoveryIsPossible
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
IECount=0
Call Get_IntersectedEdges(Dim,NC,Corn,Neib,X,Y,FrontEdges,States,Fronts,N_c,N_d,IntersectedEdges,IECount,TopEdgeRecoveryIsPossible)

if(TopEdgeRecoveryIsPossible) then

	I=0
	do
		if(I == IECount) then
			if(IECount == 0) print*,'Top Edge Recovery With NO INTERSECTION!!!'
			exit
		endif

		I = I + 1
!Part 2:		
		do J=1,NC !----------------------- Finding first Triangle having IntersectedEdge(I) as tri1 --------

			if(Corn(J,4) == 0 .And. Corn(J,1)/=Corn(J,2)) then

				do K=1,3

					if(Corn(J,K) == IntersectedEdges(I,1)) then
						if(K == 1) then
							if(Corn(J,2)==IntersectedEdges(I,2) .Or. Corn(J,3)==IntersectedEdges(I,2)) then
								tri1 = J
							endif
						elseif(K == 2) then
							if(Corn(J,1)==IntersectedEdges(I,2) .Or. Corn(J,3)==IntersectedEdges(I,2)) then
								tri1 = J
							endif
						elseif(K == 3) then
							if(Corn(J,1)==IntersectedEdges(I,2) .Or. Corn(J,2)==IntersectedEdges(I,2)) then
								tri1 = J
							endif
						endif
					endif

				end do	

			endif

		end do
!Part 3:
		do L=1,3 !----------------------- Finding second Triangle having IntersectedEdge(I) as tri2 --------

			if(Corn(tri1,L) /= IntersectedEdges(I,1) .And. Corn(tri1,L) /= IntersectedEdges(I,2)) then

				tri2 = Neib(tri1,L)
				Iface1 = L

			endif

		end do

		do M=1,3
			if( Neib(tri2,M)==tri1 ) Iface2=M 
		end do

		
		if(Iface1==1)then
		  I1=2
		  I2=3
		elseif(Iface1==2)then
		  I1=3
		  I2=1
		elseif(Iface1==3)then
		  I1=1
		  I2=2
		endif
		
		P1=Corn(tri1,I1)
		P2=Corn(tri1,I2) 
		P3=Corn(tri1,Iface1)
		P4=Corn(tri2,Iface2)
		
		if(SwapIsSafe(Dim,Corn,X,Y,tri1,tri2)) then
!Part 4:
			Call SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri1,tri2)
			print *,'swapped!!!',P3,P4
			!-------------------------- Deleting IntersectedEdges(I) from list ---------------------

			!       By incrementing loop's counter edge is deleted from the list anyway

			!-------------- Adding new intersecting edge to the end of list if exist ---------------
!Part 5:            
			if(IntersectionOccur(Dim,N_c,N_d,P3,P4,X,Y)) then
				print *,'New Edge of Swap Inserted',P3,P4
				IECount = IECount + 1
				IntersectedEdges(IECount,1) = P3
				IntersectedEdges(IECount,2) = P4
			endif
        else
!Part 6:            
			print *,'Edge Inserted to EOL',IntersectedEdges(I,1),IntersectedEdges(I,2)
			IECount = IECount + 1
			IntersectedEdges(IECount,1) = IntersectedEdges(I,1)
			IntersectedEdges(IECount,2) = IntersectedEdges(I,2)	

		endif
		
	end do

endif
print *,'Sucssesfully Exited'
!===========================================================================================
End Subroutine TopEdgeRecovery						       
!*********************************************************************************************
