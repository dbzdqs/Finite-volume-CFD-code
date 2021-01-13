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
Subroutine CalcCentroidOfQuad(Dim,Corn,X,Y,Q,x_value,y_value)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Q
Intent(Out)::x_value,y_value

Integer::Dim,Q,a,b,c,d
Integer,Dimension(1:Dim,1:4)::Corn
Real(8)::Mab_x,Mab_y,Mcd_x,Mcd_y,Mad_x,Mad_y,Mbc_x,Mbc_y,Dx0,Dy0,Dx1,Dy1,s,t,Deltax,Deltay,x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!------------- Reference: http://mathworld.wolfram.com/Quadrilateral.html ------------------
!Part 1:
a = Corn(Q,1)
b = Corn(Q,2)
c = Corn(Q,3)
d = Corn(Q,4)

!Part 2:

Mab_x = (X(a)+X(b))/2
Mab_y = (Y(a)+Y(b))/2

Mcd_x = (X(c)+X(d))/2
Mcd_y = (Y(c)+Y(d))/2

Mad_x = (X(a)+X(d))/2
Mad_y = (Y(a)+Y(d))/2

Mbc_x = (X(b)+X(c))/2
Mbc_y = (Y(b)+Y(c))/2

!Part 3:
!-------------------------------- D0: Vector From Mab to Mcd -------------------------------
Dx0 = Mcd_x - Mab_x 
Dy0 = Mcd_y - Mab_y

!-------------------------------- D1: Vector From Mad to Mbc -------------------------------
Dx1 = Mbc_x - Mad_x 
Dy1 = Mbc_y - Mad_y

!Part 4:

Deltax = Mad_x - Mab_x
Deltay = Mad_y - Mab_y

s = (Deltax*Dy1 - Dx1*Deltay)/(Dx0*Dy1 - Dx1*Dy0)
t = (Deltax*Dy0 - Dx0*Deltay)/(Dx0*Dy1 - Dx1*Dy0)

x_value = Mab_x + s*Dx0
y_value = Mab_y + s*Dy0 

!===========================================================================================
End Subroutine CalcCentroidOfQuad
!*********************************************************************************************
