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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine GeoCalAnyShape3D(Dim,NFt,IDSt,X,Y,Z,FaceTypet,Volt,Nxt,Nyt,Nzt)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFt,IDSt,X,Y,Z,FaceTypet
 Intent(Out  )::Volt,Nxt,Nyt,Nzt

 Integer::Dim,I,J,NFt,P1,P2,P3,ME,NE,NFacePnt,IFace
 Real(8)::a,b,c,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,X21,X31,Y21,Y31,Z21,Z31,DV,MagN,Volt,aa,bb,cc,x32,y32,z32
 Integer,Dimension(1:4,1:8)::IDSt
 Integer,Dimension(1:8)::FaceTypet
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:8)::Nxt,Nyt,Nzt,DAt
!*********************************************************************************************	

Volt =0.0
DAt(:) =0.0

 DO IFace=1,NFt
    NFacePnt = FaceTypet(IFace)

    DV =0.0 
    Do I=2,NFacePnt-1

       P1= IDSt(1  ,IFace)
	   P2= IDSt(I  ,IFace)
	   P3= IDSt(I+1,IFace)

	   X1 = X(P1) ; Y1 = Y(P1) ; Z1 = Z(P1)
	   X2 = X(P2) ; Y2 = Y(P2) ; Z2 = Z(P2)
	   X3 = X(P3) ; Y3 = Y(P3) ; Z3 = Z(P3)
       
       x21 = x2 - x1
       y21 = y2 - y1
       z21 = z2 - z1

       x31 = x3 - x1
       y31 = y3 - y1
       z31 = z3 - z1

       a = y21*z31 - z21 * y31
       b = z21*x31 - x21 * z31
       c = x21*y31 - y21 * x31

       MagN = Dsqrt( a*a + b*b + c*c )

	   DV = DV + a*( x1+ x2+ x3 ) + b*( y1+ y2+ y3 ) + c*( z1+ z2+ z3 )

       DAt(IFace) = DAt(IFace) + 0.5*MagN
       
    End Do 
    DV = DV /18

    NXt(IFace) = a / MagN
	NYt(IFace) = b / MagN
	NZt(IFace) = c / MagN

    NXt(IFace) = NXt(IFace) * DAt(IFace)
	NYt(IFace) = NYt(IFace) * DAt(IFace)
	NZt(IFace) = NZt(IFace) * DAt(IFace)

    Volt =  Volt + DV

 End do

!*********************************************************************************************
 End
!###########################################################################################