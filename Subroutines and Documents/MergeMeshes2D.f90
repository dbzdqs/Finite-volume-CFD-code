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
Subroutine MergeMeshes2D(Dim,NP,NF,IDS,NR,NFR,NC,NPBL,NF_BL,IDS_BL,NR_BL,NFR_BL,NC_BL,X,Y,XBL,YBL)
implicit none
!*********************************************************************************************
Intent(In   )::Dim
Intent(InOut)::NP,NF,IDS,NR,NFR,NC,NPBL,NF_BL,IDS_BL,NR_BL,NFR_BL,NC_BL,X,Y,XBL,YBL

Integer::Dim
Integer::NP
Integer::NF !Number of Faces Constructing Mesh
Integer::NC
Integer::NR
Integer::NPBL
Integer::NF_BL
Integer::NC_BL
Integer::NR_BL
Integer::FaceA,FaceB
Integer::PSA,PDA,PSB,PDB
Integer::IA,SumA,IB,JB,SumB,NMb,NFRLA,NFRLB,NIP,NIP2,PSBIsMerge,PDBIsMerge,NX,NY
Integer,Dimension(1:100)::NFR,BC,NFR_BL,BC_BL
Integer,Dimension(1:4,1:Dim)::IDS,IDS_BL
Real(8),Dimension(1:Dim)::X,Y,XBL,YBL,BoundFaceB,RegionBound,IP,IP2
Real(16)::Dif1,Dif2,Dif3,Dif4,Dif5,Dif6,Dif7,Dif8
!*********************************************************************************************
!Part 1
NMb = 1
NFRLA=1
NIP=1
BoundFaceB(:) = 0
RegionBound(:) = 0
!Part 2
Do IA = 1,NR
    Do SumA=1,NFR(IA)
        FaceA = NFRLA
        NFRLA = NFRLA + 1
        PSA = IDS(3,FaceA)
        PDA = IDS(4,FaceA)
        NFRLB = 1
        !Part 3
        Do IB=1,NR_BL
            Do SumB=1,NFR_BL(IB)
                FaceB=NFRLB
                NFRLB = NFRLB +1
                PSB=IDS_BL(3,FaceB)
                PDB=IDS_BL(4,FaceB)

                !Part 4
                Dif1 = ABS(X(PSA) - XBL(PSB));Dif2 = ABS(Y(PSA) - YBL(PSB));
                Dif3 = ABS(X(PDA) - XBL(PDB));Dif4 = ABS(Y(PDA) - YBL(PDB));
                
                Dif5 = ABS(X(PSA) - XBL(PDB));Dif6 = ABS(Y(PSA) - YBL(PDB));
                Dif7 = ABS(X(PDA) - XBL(PSB));Dif8 = ABS(Y(PDA) - YBL(PSB));
                If(Dif1+Dif2+Dif3+Dif4<0.00000000000001 .or. Dif5+Dif6+Dif7+Dif8<0.000000000000001) then
                    RegionBound(IB)=-1
                    IP(NIP) =  PDB
                    NIP=NIP+1
                    IP(NIP)=PSB
                    NIP=NIP+1
                    !Part 5
                    If(IDS(1,FaceA)==0) then
                        IDS(1,FaceA)= IDS_BL(2,FaceB)+NC
                        BoundFaceB(NMb) = FaceB
                        NMb=NMb+1
                    Else
                        IDS(2,FaceA)=IDS_BL(1,FaceB)+NC
                        BoundFaceB(NMb) = FaceB
                        NMb=NMb+1
                    End if
                End if
            End Do
            SumB=1
        End Do
        IB= 1
    End Do
    SumA =1
End Do

!Part 6
SumB=1
NFRLB=1
NIP2 = NIP
NIP=1
PSBIsMerge = 0
PDBIsMerge = 0
!Part 7
Do IB=1,NR_BL
    If(RegionBound(IB)/=-1) then
        NR = NR+1
        NFR(NR)=NFR_BL(IB)
    End if

    !Part 8
22      Do while(SumB<=NFR_BL(IB))
        JB = 1
        NIP = 1
        FaceB=NFRLB
        NFRLB=NFRLB + 1
        PSB=IDS_BL(3,FaceB)
        PDB=IDS_BL(4,FaceB)
        PSBIsMerge = 0
        PDBIsMerge = 0

        !Part 9
        Do while(BoundFaceB(JB)/=0)
            If(BoundFaceB(JB)==FaceB) then
                SumB = SumB+1
                Go to 22
            End if
            JB=JB+1
        End Do
        NF=NF+1

        !Part 10
        Do while(IP(NIP)/=0)
            If(IP(NIP)==PSB) then
                PSBIsMerge=1
            End if
            If(IP(NIP)==PDB) then
                PDBIsMerge=1
            End if
            NIP=NIP+1
        End Do

        !Part 11
        If(IDS_BL(1,FaceB)==0) then
            IDS(1,NF) = 0
        Else
            IDS(1,NF)= IDS_BL(1,FaceB)+NC
        End if
        If(IDS_BL(2,FaceB)==0) then
            IDS(2,NF)= 0
        Else
            IDS(2,NF) = IDS_BL(2,FaceB)+NC
        End if

        !Part 12
        If(PSBIsMerge==1) then
            NX=1
            Do while(1)
                Dif1 = ABS(X(NX)-XBL(PSB));Dif2 = ABS(Y(NX)-YBL(PSB));
                If(Dif1+Dif2<0.000000000000001)then
                    IDS(3,NF)=NX
                    exit
                End if
                NX = NX +1
            End Do
        Else
            NP = NP+1
            IDS(3,NF) = NP
            X(NP)=XBL(PSB)
            Y(NP)=YBL(PSB)
            IP(NIP2)=PSB
            NIP2 = NIP2+1
        End if

        !Part 13
        If(PDBIsMerge==1) then
            NX=1
            Do while(1)
                Dif1 = ABS(X(NX)-XBL(PDB));Dif2 = ABS(Y(NX)-YBL(PDB));
                If(Dif1+Dif2<0.000000000000001)then
                    IDS(4,NF) = NX
                    exit
                End if
                NX = NX +1
            End Do
        Else
            NP = NP + 1
            IDS(4,NF) = NP
            X(NP) = XBL(PDB)
            Y(NP) = YBL(PDB)
            IP(NIP2) = PDB
            NIP2 = NIP2 + 1
        End if
        SumB=SumB+1
    End Do
    SumB=1
End Do
!*********************************************************************************************
 END
!###########################################################################################