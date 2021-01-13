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
!// Developed by: *//*-+/01                        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine AddToSmoothList(Dim,Point,Elm,SmoothNode,PC)
Implicit None
!===========================================================================================
Intent(In)::Dim,Point,Elm
Intent(InOut)::SmoothNode,PC

Integer,Parameter::NODE = 1
Integer,Parameter::ELEMENT = 2

Integer::Dim,Point,Elm,PC
Integer,Dimension(1:Dim,1:2)::SmoothNode
Logical::IsInSmoothList
!===========================================================================================
!Part 1:
if(.Not. IsInSmoothList(Dim,Point,SmoothNode,PC)) then

	PC = PC + 1

	SmoothNode(PC,NODE) = Point
	SmoothNode(PC,ELEMENT) = Elm

endif
!===========================================================================================
End Subroutine AddToSmoothList
!*********************************************************************************************
