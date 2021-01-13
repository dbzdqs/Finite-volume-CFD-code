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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function isFrontEdge(Dim,Fronts,FrontEdges,States,Ei)
Implicit None
!===========================================================================================
Intent(In)::Dim,Fronts,FrontEdges,States,Ei

Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2

Integer::Dim,Fronts,I
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::FrontEdges
Logical::isFrontEdge
!===========================================================================================
!Part 1:
isFrontEdge = .false.
do I=1,Fronts
	if(States(I) /= Processed) then
		if((FrontEdges(I,1)==Ei(1) .And. FrontEdges(I,2)==Ei(2)) .Or. (FrontEdges(I,1)==Ei(2) .And. FrontEdges(I,2)==Ei(1))) then
			isFrontEdge = .true.
			exit					
		endif
	endif
end do
!===========================================================================================
End Function isFrontEdge
!*********************************************************************************************
