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
!// Date: Apr., 25, 2013                                                                   //!
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine PutPointsOnLine(Spacinge,NInitPt,NFinalPt,InitX,InitY,InitZ,XFinal,YFinal,ZFinal)
Implicit None
!*********************************************************************************************
Integer::I
Real(8)::X1,X2,Y1,Y2,Z1,Z2,LenLine,LenSegs

Integer::NInitPt
Integer::NFinalPt
Real(8),Dimension(1:NFinalPt)::InitX
Real(8),Dimension(1:NFinalPt)::InitY
Real(8),Dimension(1:NFinalPt)::InitZ
Real(8),Dimension(1:NFinalPt)::XFinal
Real(8),Dimension(1:NFinalPt)::YFinal
Real(8),Dimension(1:NFinalPt)::ZFinal
Real(8),Dimension(1:NFinalPt)::Spacinge
!*********************************************************************************************
!Part 1
X1 = InitX(1)
Y1 = InitY(1)
Z1 = InitZ(1)

X2 = InitX(2)
Y2 = InitY(2)
Z2 = InitZ(2)

XFinal(1) = X1
YFinal(1) = Y1
ZFinal(1) = Z1

LenLine = Dsqrt( (X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2 )

LenSegs = 0.0

Do I=2,NFinalPt-1

    LenSegs = LenSegs + Spacinge(I-1)

    XFinal(I) = X1 + (LenSegs/LenLine)*(X2-X1)
    YFinal(I) = Y1 + (LenSegs/LenLine)*(Y2-Y1)
    ZFinal(I) = Z1 + (LenSegs/LenLine)*(Z2-Z1)

End do

XFinal(NFinalPt) = X2
YFinal(NFinalPt) = Y2
ZFinal(NFinalPt) = Z2
!*********************************************************************************************
End
!###########################################################################################