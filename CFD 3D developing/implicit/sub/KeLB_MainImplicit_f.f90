!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate Turbulence Viscosity and Fluctuation Term of Turbulent Flow//!
!//              by Lam-Beremhorst Turbulence Model                                      //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: H. Nazari, M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                //!
!// Doc ID: MC5F101F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KeLB_MainImplicit_f(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,&
                              NP,IDS,FaceType,INW,XC,YC,ZC,DWAL,DA,DT,Vol,MR,Wb,WNP1,Mu,P,&
                              GM,DUY,InonzeroCell,InonzeroFace,NnonzeroCell,WTNP1,Mut)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DWAL,DA,DT,Vol,MR,Wb,WNP1,Mu,P,&
                &GM,DUY,InonzeroCell,InonzeroFace,NnonzeroCell
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,I,II,J,NS,ME,KK
 Real(8)::K,Epsilon,Rho,Co,tauwall,ustar,yplus,Yn,fmu
 Real(8)::Ceps1,Ceps2,SigK,Sige,CMu,Lam1,Lam2,Lam3
!-------------------------------------------------------------------------------------------
 Integer::NC
 Integer::NP
 Integer::NF
 Integer::NFW1,NFW2
 Integer::NF1,NF2
 Integer::NFI1,NFI2
 Integer::NFO1,NFO2
 Integer::NFS1,NFS2
 Integer::NFF1,NFF2
 Real(8)::MR
 Real(8)::GM,sor,DD
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW,NnonzeroCell
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTNP1v
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:Dim)::Vol
 Real(8),Dimension(1:Dim)::Mu
 Real(8),Dimension(1:Dim)::Mut
 Real(8),Dimension(1:Dim)::DT
 Real(8),Dimension(1:2,1:Dim)::WTB,WTBv
 Real(8),Dimension(1:2,1:Dim)::WTN
 Real(8),Dimension(1:2,1:Dim)::Cont,Contv
 Real(8),Dimension(1:2,1:Dim)::Dift,Diftv
 Real(8),Dimension(1:2,1:Dim)::St,Stv
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:Dim)::NX,NY,NZ
 Real(8),Dimension(1:Dim)::XC,YC,ZC
 Real(8),Dimension(1:Dim)::DA,DWAL
 Real(8),Dimension(1:Dim)::DUY
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Real(8),Dimension(1:Dim)::DKX_F,DKY_F,DKZ_F
 Real(8),Dimension(1:Dim)::DEPSX_F,DEPSY_F,DEPSZ_F
 Real(8),Dimension(1:2,1:Dim)::Rest
 Real(8),Dimension(1:2,1:Dim)::Restv

 Integer,Dimension(1:10,1:Dim)::InonzeroCell,InonzeroFace
 Real(8),Dimension(1:Dim)::eigen,JacobiL,JacobiR
 Real(8),Dimension(1:Dim)::eigenv,JacobiLv,JacobiRv
 Real(8),Dimension(1:2,1:Dim)::DW_star,DW
!*******************************************************************************************
!Part 1:
 Ceps1    =1.44D0
 Ceps2    =1.92D0
 SigK     =1.0D0
 Sige     =1.3D0
 CMu      =0.115D0
 Lam1     =1.0D0
 Lam2     =0.4D0
 Lam3     =0.2D0
 DW=0.0

!Part 2:
 Do I=1,NC
    WTN(1,I) =WTNP1(1,I)
    WTN(2,I) =WTNP1(2,I)
 End Do


!Part 3:
Call KeLam_BC3D(Dim,NX,NY,NZ,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB)

!Part 4:
Call Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)

!Part 5:
Call KFi_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1&
                    &,WTNP1,WB,WTB,DKX_F,DKY_F,DKZ_F,DEPSX_F,DEPSY_F,DEPSZ_F)

!Part 6:
Call KFi_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)

!Part 7:
Call KFi_Dif3d(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DEPSX_F,DEPSY_F,DEPSZ_F,Sigk,Sige,MR,Mu,Mut,Dift)

!Part 8:
Call KeLB_Source3D(Dim,NC,IDS,DWAL,Vol,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,&
                 GM,P,WTNP1,Mu,Mut,DUY,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,St)

 Rest(1,:) = Cont(1,:) + Dift(1,:) + St(1,:)
 Rest(2,:) = Cont(2,:) + Dift(2,:) + St(2,:)

!Part 9.1:
Call Calculate_eigValTurb(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA,eigen)

!Part 9.2:
Call IncrementTurb(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA,JacobiL,JacobiR)

sor=0.2
!!!===============================================================================================
!!Part 10:
Call LUSGSTurb(Dim,NC,NX,NY,NZ,DA,GM,Rest,NnonzeroCell,InonzeroCell,InonzeroFace,eigen,JacobiL,JacobiR,Vol,DT,WTNP1,DW)

!Part 11:
Do J=1,NC

  !part 12:
   WTNP1(1,J) = WTN(1,J) + sor*DW(1,j)
   WTNP1(2,J) = WTN(2,J) + sor*DW(2,j)

  !Part 13:
   if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
   if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)

  !Part 14:
   Rho     = WNP1(1,J)
   K       = WTNP1(1,J)/Rho
   Epsilon = WTNP1(2,J)/Rho

  !Part 15:
   II = INW(J)     ! II: Wall Face
   Yn = DWAL(J)
   ME = IDS(1,II)

   Tauwall = Mu(ME)*DUY(II) !!!!!!!!!!!!!!!!!!!!
   Ustar   = Dsqrt(abs(MR*Tauwall/WNP1(1,ME)))
   Yplus   = (1.0/MR)*Rho*Ustar*Yn/Mu(J)

   FMu = 0.04 + (1.0 - 0.04)*((1.0 - EXP(-(Yplus - 8.0)/26.0))**2)

  !Part 16:
   Mut(J) =(1.0/MR)*Cmu*FMU*Rho*K*K / (Epsilon+1.e-20)

End Do !J

!###########################################################################################
 End
!*******************************************************************************************
