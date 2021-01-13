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
Subroutine CheckContinuity(Dim,Fronts,FrontEdges,States,Num_Of_Loops,Continuity)
Implicit None
!===========================================================================================
Intent(In)::Dim,Fronts,FrontEdges,States
Intent(Out)::Num_Of_Loops,Continuity

Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2

Integer::Dim,Fronts,ActiveFronts,Num_Of_Loops,Start,Next,index,I,J,C,FC,cnt
Integer,Dimension(1:Dim)::States,Loops
Integer,Dimension(1:Dim,1:4)::FrontEdges
logical::AlreadyProcessed,Continuity,LoopCompleted
!===========================================================================================
Continuity = .True.
ActiveFronts = 0
Num_Of_Loops = 0
FC = 0
!Part 1:
do I=1,Fronts
	if(States(I) /= Processed) then
		ActiveFronts  = ActiveFronts + 1
	endif
end do

if(ActiveFronts > 0) then

!Part 2:

do
	
	do I=1,Fronts

		if(States(I) /= Processed) then

			AlreadyProcessed = .False.

			do J=1,FC
				if(Loops(J) == I) then
					AlreadyProcessed = .True.
					exit
				endif
			end do

			if(.Not. AlreadyProcessed) then

				Start = FrontEdges(I,LeftVertex)
				index = I
				exit

			endif 

		endif

	end do

	Next = Start
	LoopCompleted = .False.
	C = 0
	cnt = 0
!Part 3:
	do
		do I=1,Fronts
			if(States(I) /= Processed) then

				AlreadyProcessed = .False.

				do J=1,FC
					if(Loops(J) == I) then
						AlreadyProcessed = .True.
						exit
					endif
				end do

				if(.Not. AlreadyProcessed) then
					
					if(FrontEdges(I,RightVertex) == Next) then
						
						FC = FC + 1
						Loops(FC) = I
						cnt = cnt + 1
						Next = FrontEdges(I,LeftVertex)
						
						if(Next == Start) then
							Num_Of_Loops = Num_Of_Loops + 1
							LoopCompleted = .True.
							print *,'Loop: ',Num_Of_Loops,' : ',cnt
						endif
						exit
					endif

				endif

			endif

		end do

		C = C + 1

		if(C > ActiveFronts) exit

		if(LoopCompleted) exit

	end do
!Part 4:
	if(C > ActiveFronts) then
		Continuity = .False.
		print *,'No Match for: ',Next
		print *,'>>>>>>>>>>>>>>> Loop Is BROKEN!!! <<<<<<<<<<<<<<<<<'
		exit
	endif

	if(FC == ActiveFronts) exit

end do

else
	print *,'>>>>>>>>>>>>>>> Done!!! <<<<<<<<<<<<<<<'
endif
!===========================================================================================
End Subroutine CheckContinuity
!*********************************************************************************************
