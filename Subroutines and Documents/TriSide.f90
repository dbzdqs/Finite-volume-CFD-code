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
Function TriSide(Dim,Corn,X,Y,Tri)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Tri

Integer::Dim,Tri,A,B,C,TriSide
Integer,Dimension(1:Dim,1:4)::Corn
Real(8)::Ux,Uy,Vx,Vy,CrossProduct
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
A = Corn(Tri,1)
B = Corn(Tri,2)
C = Corn(Tri,3)

Ux = X(C) - X(A)
Uy = Y(C) - Y(A)
Vx = X(C) - X(B)
Vy = Y(C) - Y(B)

CrossProduct = Ux*Vy - Uy*Vx

TriSide = INT(CrossProduct/DABS(CrossProduct))

!===========================================================================================
End Function TriSide
!*********************************************************************************************
