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
Subroutine ExtrudeLineMesh(Dim,NCurv,NedgCurvs,NBL,NEdgCurv,NPtCurvs,EdgPt,IConectedEdg,BLPt,CuspPt,NC_BL,NF_BL,NR_BL,NFR_BL,IDS_BL)
Implicit None
!*********************************************************************************************
Intent(In   )::Dim,NCurv,NedgCurvs,NBL,NEdgCurv,NPtCurvs,EdgPt,IConectedEdg,CuspPt
Intent(Out  )::NC_BL,NF_BL,NR_BL,NFR_BL,IDS_BL
Intent(InOut)::BLPt

Integer::Dim,I,J,J1,JJ,NC_BL,NF_BL,NR_BL,Sum,E1,E2,Pt1,Pt2,P1,P2,ME,NE
Integer,Dimension(1:100)::NFR_BL
Integer,Dimension(1:4,1:Dim)::IDS_BL
Integer::NCurv
Integer::NedgCurvs
Integer::NPtCurvs
Integer::NBL
Integer,Dimension(1:100)::NEdgCurv
Integer,Dimension(1:2,1:Dim)::EdgPt
Integer,Dimension(1:2,1:Dim)::IConectedEdg
Integer,Dimension(1:25,1:Dim)::BLPt
Integer,Dimension(1:Dim)::CuspPt
!*********************************************************************************************
!Part 1
NR_BL = 0
NF_BL = 0
Sum = 0
Do I=1,NCurv
    NR_BL = NR_BL + 1
    NFR_BL(NR_BL) = NEdgCurv(I)
    
    !Part 2
    Do J=Sum+1,Sum+NEdgCurv(I)
        Pt1 = EdgPt(1,J)
        Pt2 = EdgPt(2,J)

        P1 = BLPt(1,Pt1)
        P2 = BLPt(1,Pt2)

        ME = J
        NE = 0
        
        NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2
    End Do
    Sum=Sum+NEdgCurv(I)

End Do
    
!Part 3
Do J=1,NPtCurvs
    IF( CuspPt(J)==1 )Then
        Do I=2,NBL+1
            BLPt(I,J) = BLPt(1,J)
        End Do
        Cycle
    EndIF
    E1 = IConectedEdg(1,J)
    E2 = IConectedEdg(2,J)
        
    !Part 4
    IF(E1==0)Then
        NR_BL = NR_BL + 1
        NFR_BL(NR_BL) = NBL
        Do J1=1,NBL
            P1 = BLPt(J1+1,J)
            P2 = BLPt(J1  ,J)

            ME = E1 + (J1-1)*NEdgCurvs+1
            NE = 0

            NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2
        End Do
    !Part 5
    ElseIF(E2==0)Then
        NR_BL = NR_BL + 1
        NFR_BL(NR_BL) = NBL
        Do J1=1,NBL
            P1 = BLPt(J1  ,J)
            P2 = BLPt(J1+1,J)

            ME = E2 + J1*NEdgCurvs
            NE = 0

            NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2
        End Do
    EndIF
End Do
    
!Part 6
NR_BL = NR_BL + 1
Do J=1,NPtCurvs
    IF( CuspPt(J)==1 )Then
        Do I=2,NBL+1
            BLPt(I,J) = BLPt(1,J)
        End Do
        Cycle
    EndIF

    E1 = IConectedEdg(1,J)
    E2 = IConectedEdg(2,J)
        
    !Part 7
    IF(E1/=0 .and. E2/=0)Then
        NFR_BL(NR_BL) = NFR_BL(NR_BL) + NBL
        Do J1=1,NBL
            P1 = BLPt(J1  ,J)
            P2 = BLPt(J1+1,J)

            ME = E1 + (J1-1)*NEdgCurvs
            NE = E2 + (J1-1)*NEdgCurvs

            NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2
        End Do
    EndIF
End Do
    
!Part 8
NFR_BL(NR_BL) = NFR_BL(NR_BL) + (NBL-1)*NEdgCurvs
Do I=2,NBL
    Do J=1,NEdgCurvs

        Pt1 = EdgPt(1,J)
        Pt2 = EdgPt(2,J)

        P1 = BLPt(I,Pt1)
        P2 = BLPt(I,Pt2)

        ME = J + (I-1) * NEdgCurvs
        NE = J + (I-2)   *NEdgCurvs

        NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2

    End Do
End Do
    
!Part 9
Sum = 0
Do I=1,NCurv

    NR_BL = NR_BL + 1
    NFR_BL(NR_BL) = NEdgCurv(I)

    Do J=Sum+1,Sum+NEdgCurv(I)

        Pt1 = EdgPt(1,J)
        Pt2 = EdgPt(2,J)

        P2 = BLPt(NBL+1,Pt1)
        P1 = BLPt(NBL+1,Pt2)

        ME = J+ (NBL-1) * NEdgCurvs
        NE = 0

        NF_BL = NF_BL+1 ; IDS_BL(1,NF_BL) = ME ; IDS_BL(2,NF_BL) = NE ; IDS_BL(3,NF_BL) = P1 ; IDS_BL(4,NF_BL) = P2

    End Do

    Sum=Sum+NEdgCurv(I)
End Do

!Part 10
 NC_BL = NedgCurvs*NBL
!*********************************************************************************************
END
!##########################################################################################