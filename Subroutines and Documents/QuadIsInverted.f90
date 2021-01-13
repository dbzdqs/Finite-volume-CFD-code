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
Function QuadIsInverted(Dim,Corn,X,Y,Quad)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Quad

Integer::Dim,Quad,C,I
Integer,Dimension(1:4)::QSide
Integer,Dimension(1:Dim,1:4)::Corn 
Logical::QuadIsInverted
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
QuadIsInverted = .False.

C = 0
!Part 1:
Call QuadSide(Dim,Corn,X,Y,Quad,QSide)
	
do I=1,4
	if(QSide(I)==-1) C=C+1
end do

if(C>1) then

	QuadIsInverted = .True.

endif

!===========================================================================================
End Function QuadIsInverted
!*********************************************************************************************
