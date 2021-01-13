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
Subroutine CalcBisectorIntersection(Dim,X,Y,A,B,C,x_value,y_value)
Implicit None
!===========================================================================================
Intent(In)::Dim,X,Y,A,B,C
Intent(Out)::x_value,y_value

Integer::Dim,A,B,C
Real(8)::Ux,Uy,Vx,Vy,norm1,norm2,BVx,BVy,value,t,Deltax,Deltay,x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--- Notice: this subroutine calcultates the coordinate of intersection point of Bisector Vector of AB and AC with BC
!Part 1:
Ux = X(B) - X(A) !--------------- Defining Vector U (A to B) ------------------
Uy = Y(B) - Y(A)
            
Vx = X(C) - X(A) !-------------- Defining Vector V (A to C) -------------------
Vy = Y(C) - Y(A)

!Part 2:

norm1 = DSQRT(Ux*Ux + Uy*Uy) !------------- Norm of U vector -------------------
norm2 = DSQRT(Vx*Vx + Vy*Vy) !------------- Norm of V vector -------------------

!Part 3:

BVx = norm1*Vx + norm2*Ux !------------ Defining Bisector Vector ------------- 
BVy = norm1*Vy + norm2*Uy

Ux = X(C) - X(B) !------------ Defining Vector U (from B to C) --------------
Uy = Y(C) - Y(B) 
               
!--- Checking Intersection of Bisector Vector (BV) and current Vector (U) ---
!----------- Consider BV = D0 as a vector in form of (P0 + sD0) -------------
!----------- Consider U = D1 as a vector in form of (P1 + tD1) --------------
!----------------- According to that P0 = A and P1 = B ---------------------

!Part 4:

Deltax = X(B) - X(A)  
Deltay = Y(B) - Y(A)
               
value = BVx*Uy - BVy*Ux
               
t = (Deltax*BVy - Deltay*BVx)/value

x_value = X(B) + t*Ux 
y_value = Y(B) + t*Uy
!===========================================================================================
End Subroutine CalcBisectorIntersection 
!*********************************************************************************************
