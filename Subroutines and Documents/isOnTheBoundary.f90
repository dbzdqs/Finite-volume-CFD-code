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
Function isOnTheBoundary(Dim,Fronts,point,FrontEdges,States)
Implicit None
!===========================================================================================
Intent(In)::Dim,Fronts,point,FrontEdges,States

Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2

Integer::Dim,Fronts,point,I
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::FrontEdges
Logical::isOnTheBoundary
!===========================================================================================
!Part 1:
isOnTheBoundary = .False.
do I=1,Fronts
	if(States(I) /= -1) then
		if(FrontEdges(I,LeftVertex) == point .Or. FrontEdges(I,RightVertex) == point) then
			isOnTheBoundary = .True.
			exit	
		endif
	endif
end do
!===========================================================================================
End Function isOnTheBoundary
!*********************************************************************************************
