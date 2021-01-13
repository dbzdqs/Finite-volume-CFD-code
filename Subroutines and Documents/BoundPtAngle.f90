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
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine BoundPtAngle(Dim,NPtCurvs,IConectedEdg,EdgPt,XBL,YBL,AngleBP)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,NPtCurvs,IConectedEdg,EdgPt,XBL,YBL
Intent(Out   )::AngleBP

Integer::Dim,I,J,P1,P2,P3,Pt,E1,E2
Real(8)::TOP,DoWN,PI,RAD_TO_DEG
Real(8),Dimension(1:2)::L1,L2
Integer::NPtCurvs
Integer,Dimension(1:2,1:Dim)::EdgPt
Integer,Dimension(1:2,1:Dim)::IConectedEdg
Real(8),Dimension(1:Dim)::XBL,YBL
Real(8),Dimension(1:Dim)::AngleBP
!*********************************************************************************************
!Part 1
PI = 4*Atan(1.0)
RAD_TO_DEG  = 	180/PI

!Part 2
Do J=1,NPtCurvs
    E1 = IConectedEdg(1,J)
    E2 = IConectedEdg(2,J)

    !Part 3
    IF(E1==0 .or. E2==0)Then
        AngleBP(J) = 0.0
    Else
        !Part 4
        P1 = EdgPt(1,E1)
        P2 = EdgPt(2,E1)
        P3 = EdgPt(2,E2)

        L1(1) = XBL(p2)-XBL(p1) ; L1(2) = YBL(p2)-YBL(p1)
        L2(1) = XBL(p3)-XBL(p1) ; L2(2) = YBL(p3)-YBL(p1)

        TOP = L1(1)*L2(1)+L1(2)*L2(2)
        DoWN= DSQRT(L1(1)**2+L1(2)**2)*DSQRT(L2(1)**2+L2(2)**2)

        !Part 5
        IF( dabs(Down)<0.00000000001D0)  Down=0.000000000001D0

        AngleBP(J)=DACOS(TOP/DoWN)*RAD_TO_DEG

    EndIF
End Do
!*********************************************************************************************
END
!###########################################################################################