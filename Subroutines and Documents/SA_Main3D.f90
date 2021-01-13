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
!// Date: Oct., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine SA_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                      NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,WB,WNP1,&
                      Mu0,Mu,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                      NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,WB,WNP1,&
                      Mu0,Mu,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0
 Intent(InOut)::WTNP1,Mut

 Integer::Dim,I,J,Ii,Jj,NC,NF,NFW1,NFW2,NF1,NF2,NS,Ncyc,NP,NRKS,ME,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::Mu0,RKco,MR,Co
 Real(8)::Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0,Cv13,Chi,Chi3,Fv1
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,WTN,Cont,Prod,Dift,Dest
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Real(8),Dimension(1:Dim)::DNuX,DNuY,DNuZ,DRNuX,DRNuY,DRNuZ,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc
 Real(8),Dimension(1:Dim)::X,Y,Z,NX,NY,NZ,Vol,Mu,Mut,XC,YC,ZC,DW,DT
 Real(8),Dimension(1:5)::RKJ
 Integer,Dimension(1:Dim)::FaceType
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    WTN(1,I) = WTNP1(1,I)
 End Do

!Part 2:
 Do NS=1,NRKS

   !Part 3:
    Call SA_BC3D(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,NX,NY,NZ,WB,WTNP1,Nuhat0,WTB)

   !Part 4:
    Call Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
    
    Call SA_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc)

   !Part 5:
    Call SA_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,NX,NY,NZ,&
                       WNP1,WTNP1,WB,WTB,DNuX,DNuY,DNuZ)

   !Part 6:
    Call SA_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)

   !Part 7:
    Call SA_Dif3D(Dim,NC,NF1,NF2,NF,MR,Cb2,Sigma,NX,NY,NZ,Vol,WTNP1,WNP1,Mu,WTB,IDS,&
                  DNuX,DNuY,DNuZ,DRNuX,DRNuY,DRNuZ,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc,Dift)

   !Part 8:
    Call SA_ProdDest3D(Dim,NC,Vol,Dw,MR,Cv1,Cb1,Cw1,Cw2,Cw3,Kei,WNP1,Mu,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,DWY_C,DWZ_C,WTNP1,Dest,Prod)

   !Part 9:
    RKco=RKJ(NS)
    DO I=1,NC

	   Co = RKco*DT(I)/Vol(I)
       WTNP1(1,I) = WTN(1,I) - Co*( Cont(1,I) - Prod(1,I) - Dift(1,I) + Dest(1,I))

       IF( WTNP1(1,I)<0. ) WTNP1(1,I) = WTN(1,I)
    End do

   !Part 10:
    Cv13=Cv1*Cv1*Cv1
    DO I=1,NC

	   Chi  = WTNP1(1,I)/Mu(I)
	   Chi3 = Chi*Chi*Chi

	   Fv1  = Chi3/(Chi3+Cv13)

	   Mut(I) = Fv1*WTNP1(1,I)
    End do

 End do
!**********************************************************************************************
 End
!##############################################################################################






