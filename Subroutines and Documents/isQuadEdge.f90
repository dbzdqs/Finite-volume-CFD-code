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
Function isQuadEdge(p1,p2,newQuad)
Implicit None
!===========================================================================================
Intent(In)::p1,p2,newQuad

Integer::p1,p2
Integer,Dimension(1:4)::newQuad
Logical::isQuadEdge
!===========================================================================================
!Part 1:
if((newQuad(1)==p1 .And. newQuad(2)==p2) .or. (newQuad(1)==p2 .And. newQuad(2)==p1)) then
	isQuadEdge = .true.
elseif((newQuad(2)==p1 .And. newQuad(3)==p2) .or. (newQuad(2)==p2 .And. newQuad(3)==p1)) then
	isQuadEdge = .true.
elseif((newQuad(3)==p1 .And. newQuad(4)==p2) .or. (newQuad(3)==p2 .And. newQuad(4)==p1)) then
	isQuadEdge = .true.
elseif((newQuad(4)==p1 .And. newQuad(1)==p2) .or. (newQuad(4)==p2 .And. newQuad(1)==p1)) then
	isQuadEdge = .true.
else
	isQuadEdge = .false.
endif
!===========================================================================================
End Function isQuadEdge
!*********************************************************************************************
