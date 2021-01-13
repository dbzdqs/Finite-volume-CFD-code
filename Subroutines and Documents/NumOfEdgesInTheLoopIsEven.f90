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
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LStart,LEnd,Situation)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LStart,LEnd,Situation

Integer,Parameter::Processed=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,Fronts,FrontIndex,LStart,LEnd,Situation,C,I,J,Index,Prev,Next,P,Q,TEC
Integer,Dimension(1:1000)::TElms
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::NumOfEdgesInTheLoopIsEven,Match
!===========================================================================================

Select Case(Situation)

	Case(1)

		C = 0
		Q = LStart
		Prev = FrontEdges(FrontIndex,Element)

		do I=1,2
			if(FrontEdges(FrontIndex,I) /= LStart) then
				P = FrontEdges(FrontIndex,I)
				exit 
			endif
		end do

		do
			TEC = 0

			do
				do J=1,3
					if(Corn(Prev,J) == P) then
						Next = Neib(Prev,J)
						exit
					endif
				end do
				
				TEC = TEC + 1
				TElms(TEC) = Prev

				if(Next == 0) then
                    exit    
                elseif(Corn(Next,4) /= 0) then
                    exit    
                endif
				
				do J=1,3
					if(Corn(Prev,J)/=Q .And. Corn(Prev,J)/=P) then
						P = Corn(Prev,J)
						exit	
					endif
				end do

				Prev = Next
						
			end do

			Match = .False.

			do I=1,Fronts
				if(States(I) /= Processed) then
					do J=1,TEC
						if(FrontEdges(I,RightVertex) == Q .And. FrontEdges(I,Element) == TElms(J)) then		
							Match = .True.					
							Index = I
							exit	
						endif
					end do

					if(Match) exit
				endif
			end do

			P = Q
			Q = FrontEdges(Index,LeftVertex)
			Prev = FrontEdges(Index,Element) 
			C = C + 1

			if(Q == LEnd .Or. Q == LStart) exit

		end do

		if(Q == LStart) then ! ------ Case When Nm is not On the Current Loop of Fronts ------
		
			NumOfEdgesInTheLoopIsEven = .True.
		
		else ! ------------------ Case When Nm is On the Current Loop of Fronts -------------------

			if(MOD(C+1,2) == 0) then
				NumOfEdgesInTheLoopIsEven = .True.
			else
				NumOfEdgesInTheLoopIsEven = .False.
			endif

		endif

	Case(2)

		C = 0
		Q = LStart
		Prev = FrontEdges(FrontIndex,Element)

		do I=1,2
			if(FrontEdges(FrontIndex,I) /= LStart) then
				P = FrontEdges(FrontIndex,I)
				exit 
			endif
		end do

		do
			TEC = 0

			do
				do J=1,3
					if(Corn(Prev,J) == P) then
						Next = Neib(Prev,J)
						exit
					endif
				end do

				TEC = TEC + 1
				TElms(TEC) = Prev
				
				if(Next == 0) then
                    exit    
                elseif(Corn(Next,4) /= 0) then
                    exit    
                endif
				
				do J=1,3
					if(Corn(Prev,J)/=Q .And. Corn(Prev,J)/=P) then
						P = Corn(Prev,J)
						exit	
					endif
				end do

				Prev = Next
						
			end do

			Match = .False.

			do I=1,Fronts
				if(States(I) /= Processed) then
					do J=1,TEC
						if(FrontEdges(I,LeftVertex) == Q .And. FrontEdges(I,Element) == TElms(J)) then
							Index = I
							Match = .True.
						endif
					end do

					if(Match) exit

				endif
			end do

			P = Q
			Q = FrontEdges(Index,RightVertex)
			Prev = FrontEdges(Index,Element)
			C = C + 1

			if(Q == LEnd .Or. Q == LStart) exit

		end do

		if(Q == LStart) then ! ------ Case When Nm is not On the Current Loop of Fronts ------
		
			NumOfEdgesInTheLoopIsEven = .True.
			
		else ! ------------------ Case When Nm is On the Current Loop of Fronts -------------------

			if(MOD(C+1,2) == 0) then
				NumOfEdgesInTheLoopIsEven = .True.
			else
				NumOfEdgesInTheLoopIsEven = .False.
			endif

		endif
		  
End Select  
!===========================================================================================
End Function NumOfEdgesInTheLoopIsEven
!*********************************************************************************************
