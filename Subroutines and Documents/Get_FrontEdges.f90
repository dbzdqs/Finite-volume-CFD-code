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
Subroutine Get_FrontEdges(Dim,Corn,Neib,NC,FrontEdges,Fronts,The_Level) ! Tested
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,NC
Intent(Out)::FrontEdges,Fronts
Intent(Inout)::The_Level

Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4

Integer::Dim,I,J,NC,Count,Fronts,The_Level
Integer,Dimension(1:Dim,1:4)::FrontEdges,Corn,Neib
Integer,Dimension(1:3)::TriNeibs
!===========================================================================================
do I=1,NC

	Count=0	                                           ! A Counter to calculate number of triangular neibours of an element
	                                 
	if(Corn(I,4)==0) then                              ! Checking an element to be a triangular
!Part 1:		
		do J=1,3                                       ! Since triangular elements only got three neibours
			if(Neib(I,J) /= 0) then 				   ! If current neibour exists
				if(Corn(Neib(I,J),4)==0) then          ! Checking triangularity of neibours
					Count=Count+1
					TriNeibs(Count)=Neib(I,J)          ! Here we save triangle neibours
				endif
			endif
		end do

		if(Count<=2) then
!Part 2:
			do J=1,3
				
				if(Neib(I,J) == 0) then                !If opposite edge to current node is a boundary edge
					
					Fronts=Fronts+1
!Part 3:
					if(J == 1) then

						FrontEdges(Fronts,LeftVertex)=Corn(I,2)  ! Saving the LEFT corners of the discovered front edge
						FrontEdges(Fronts,RightVertex)=Corn(I,3) ! Saving the RIGHT corners of the discovered front edge

					elseif(J == 2) then

						FrontEdges(Fronts,LeftVertex)=Corn(I,3)  ! Saving the LEFT corners of the discovered front edge
						FrontEdges(Fronts,RightVertex)=Corn(I,1) ! Saving the RIGHT corners of the discovered front edge

					elseif(J == 3) then

						FrontEdges(Fronts,LeftVertex)=Corn(I,1)  ! Saving the LEFT corners of the discovered front edge
						FrontEdges(Fronts,RightVertex)=Corn(I,2) ! Saving the RIGHT corners of the discovered front edge

					endif
!Part 4:
					FrontEdges(Fronts,Element)=I           ! Saving cell number of triangle elements having a Front Edge
					FrontEdges(Fronts,Level) = The_Level   ! Saving Level of current Front Edge

				endif

			end do

		endif

	endif

end do
!===========================================================================================
End Subroutine Get_FrontEdges 
!*********************************************************************************************
