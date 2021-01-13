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
Subroutine InvSlopOfEdges(Dim,XBL,YBL,NedgCurvs,NPtCurvs,IConectedEdg,EdgPt,EdgInvSlop)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,XBL,YBL,NedgCurvs,NPtCurvs,IConectedEdg,EdgPt
Intent(Out  )::EdgInvSlop

Integer::Dim,I,J,P1,P2,Pt,E1,E2
Integer::NedgCurvs
Integer::NPtCurvs
Real(8)::SlopX,SlopY,Mag
Real(8),Dimension(1:Dim)::XBL,YBL
Real(8),Dimension(1:2,1:Dim)::EdgInvSlop
Integer,Dimension(1:2,1:Dim)::EdgPt
Integer,Dimension(1:2,1:Dim)::IConectedEdg
!*********************************************************************************************
EdgInvSlop(:,:) =  0.0

!Part 1
 Do I=1,NedgCurvs

    P1 = EdgPt(1,I)
    P2 = EdgPt(2,I)

    SlopX =-( YBL(P2) - YBL(P1) )
    SlopY = ( XBL(P2) - XBL(P1) )

    EdgInvSlop(1,P1) = EdgInvSlop(1,P1) + SlopX
    EdgInvSlop(2,P1) = EdgInvSlop(2,P1) + SlopY

    EdgInvSlop(1,P2) = EdgInvSlop(1,P2) + SlopX
    EdgInvSlop(2,P2) = EdgInvSlop(2,P2) + SlopY

 End Do

!Part 2
 Do I=1,NPtCurvs

    E1 = IConectedEdg(1,I)
    E2 = IConectedEdg(2,I)

    IF( E1/=0 .and. E2/=0 )Then
        EdgInvSlop(1,I) = EdgInvSlop(1,I) * 0.5
        EdgInvSlop(2,I) = EdgInvSlop(2,I) * 0.5
    EndIF

    Mag = Dsqrt( EdgInvSlop(1,I)**2 + EdgInvSlop(2,I)**2 )

    EdgInvSlop(1,I) = EdgInvSlop(1,I) / Mag
    EdgInvSlop(2,I) = EdgInvSlop(2,I) / Mag

 End Do
!*********************************************************************************************
 END
!###########################################################################################