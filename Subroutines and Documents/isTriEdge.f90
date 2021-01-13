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
Function isTriEdge(Q1,Q2,P1,P2,P3)
Implicit None
!===========================================================================================
Intent(In)::Q1,Q2,P1,P2,P3

Integer::P1,P2,P3,Q1,Q2
Logical::isTriEdge
!===========================================================================================
!Part 1:
if((P1==Q1 .And. P2==Q2) .or. (P1==Q2 .And. P2==Q1)) then
	isTriEdge = .true.
elseif((P2==Q1 .And. P3==Q2) .or. (P2==Q2 .And. P3==Q1)) then
	isTriEdge = .true.
elseif((P3==Q1 .And. P1==Q2) .or. (P3==Q2 .And. P1==Q1)) then
	isTriEdge = .true.
else
	isTriEdge = .false.
endif
!===========================================================================================
End Function isTriEdge
!*********************************************************************************************
