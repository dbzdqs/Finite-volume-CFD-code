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
!// Supervisor:                                                                            //!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Date: April, 01, 2017                                                                  //!
!// Supervisor:                                                                            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine GenBLayerPt2D(Dim,NPtCurv,NCurv,DistributionType,BLThick,NBL,BLPt,NPBL,XBL,YBL)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,NPtCurv,NCurv,DistributionType,BLThick,NBL
Intent(InOut)::BLPt,NPBL,XBL,YBL

Integer::Dim,I,J,JJ,P1,P2,Pt,P,Sum
Integer::NPtCurvs
Integer::NCurv
Integer::NBL
Integer::NPBL
Integer::Flag
Integer,Dimension(1:100)::NPtCurv
Integer,Dimension(100)::DistributionType
Integer,Dimension(1:25,1:Dim)::BLPt
Real(8),Dimension(1:Dim)::XBL,YBL
Real(8),Dimension(1:Dim)::BLThick
Real(8)::Lambda, LenOvr

Real(8),Dimension(1:NBL+1)::InitX
Real(8),Dimension(1:NBL+1)::InitY
Real(8),Dimension(1:NBL+1)::InitZ
Real(8),Dimension(1:NBL+1)::XFinal
Real(8),Dimension(1:NBL+1)::YFinal
Real(8),Dimension(1:NBL+1)::ZFinal
Real(8),Dimension(1:NBL+1)::Spacinge,UnitSpacinge
!*********************************************************************************************
!Part 1:
 Sum = 0
 Do J=1,NCurv

   !Part 2:
    Lambda = 0.2
    Flag = DistributionType(J)
    Call DistributionFunction(Flag,NBL,Lambda,UnitSpacinge)
        
   !Part 3:
    Do I=Sum+1,Sum+NPtCurv(J)
        
      !Part 3:
       P1 = BLPt(1    ,I)
       P2 = BLPt(NBL+1,I)

      !Part 4:
       Do JJ=1,NBL+1
          Spacinge(JJ) = BLThick(P1) * UnitSpacinge(JJ)
       End Do
             
      !Part 5:
       InitX(1) = XBL(P1)
       InitY(1) = YBL(P1)
       InitZ(1) = 0.0

       InitX(2) = XBL(P2)
       InitY(2) = YBL(P2)
       InitZ(2) = 0.0
        
      !Part 6:
       Call PutPointsOnLine(Spacinge,2,NBL+1,InitX,InitY,InitZ,XFinal,YFinal,ZFinal)

      !Part 7:
       Do JJ=2,NBL
          NPBL=NPBL+1
          XBL(NPBL) = XFinal(JJ)
          YBL(NPBL) = YFinal(JJ)
          BLPt(JJ,I) = NPBL
       End do  
    
    End Do
    Sum=Sum+NPtCurv(J)
 End Do
!*********************************************************************************************
 End
!###########################################################################################