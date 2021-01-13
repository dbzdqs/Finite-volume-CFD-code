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
Subroutine UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,tri)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Fronts,States,tri
Intent(Inout)::FrontEdges

Integer,Parameter::Processed=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,Fronts,tri,I,J
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:2)::E
Integer,Dimension(1:Dim,1:4)::Corn,FrontEdges
!===========================================================================================
!Part 1:
do I=1,3
	if(I/=3) then
		E(1) = Corn(tri,I)
		E(2) = Corn(tri,I+1)
		do J=1,Fronts
			if(States(J)/=Processed) then
				if((FrontEdges(J,LeftVertex)==E(1) .And. FrontEdges(J,RightVertex)==E(2)) .Or. (FrontEdges(J,LeftVertex)==E(2) .And. FrontEdges(J,RightVertex)==E(1))) then
					FrontEdges(J,Element) = tri
				endif
			endif
		end do
	else
		E(1) = Corn(tri,I)
		E(2) = Corn(tri,1)
		do J=1,Fronts
			if(States(J)/=Processed) then
				if((FrontEdges(J,LeftVertex)==E(1) .And. FrontEdges(J,RightVertex)==E(2)) .Or. (FrontEdges(J,LeftVertex)==E(2) .And. FrontEdges(J,RightVertex)==E(1))) then
					FrontEdges(J,Element) = tri
				endif
			endif
		end do
	endif
end do
!===========================================================================================
End Subroutine UpdateFrontElements
!*********************************************************************************************
