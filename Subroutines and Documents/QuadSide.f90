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
Subroutine QuadSide(Dim,Corn,X,Y,Quad,QSide)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Quad
Intent(Out)::QSide

Integer::Dim,Quad,A,B,C,D
Integer,Dimension(1:4)::QSide
Integer,Dimension(1:Dim,1:4)::Corn
Real(8)::Ux,Uy,Vx,Vy,CrossProduct
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
A = Corn(Quad,1)
B = Corn(Quad,2)
C = Corn(Quad,3)
D = Corn(Quad,4)
!------------------------------ Triangles ABC and ABD --------------------------------------
Ux = X(C) - X(A)
Uy = Y(C) - Y(A)
Vx = X(C) - X(B)
Vy = Y(C) - Y(B)
CrossProduct = Ux*Vy - Uy*Vx
if(CrossProduct /= 0) then
	QSide(1) = INT(CrossProduct/DABS(CrossProduct))
else
	QSide(1) = 1
endif		

Ux = X(D) - X(A)
Uy = Y(D) - Y(A)
Vx = X(D) - X(B)
Vy = Y(D) - Y(B)
CrossProduct = Ux*Vy - Uy*Vx
if(CrossProduct /= 0) then
	QSide(2) = INT(CrossProduct/DABS(CrossProduct))
else
	QSide(2) = 1
endif
!------------------------------ Triangles CDA and CDB --------------------------------------
Ux = X(A) - X(C)
Uy = Y(A) - Y(C)
Vx = X(A) - X(D)
Vy = Y(A) - Y(D)
CrossProduct = Ux*Vy - Uy*Vx
if(CrossProduct /= 0) then
	QSide(3) = INT(CrossProduct/DABS(CrossProduct))
else
	QSide(3) = 1
endif

Ux = X(B) - X(C)
Uy = Y(B) - Y(C)
Vx = X(B) - X(D)
Vy = Y(B) - Y(D)
CrossProduct = Ux*Vy - Uy*Vx
if(CrossProduct /= 0) then
	QSide(4) = INT(CrossProduct/DABS(CrossProduct))
else
	QSide(4) = 1
endif

!===========================================================================================
End Subroutine QuadSide
!*********************************************************************************************
