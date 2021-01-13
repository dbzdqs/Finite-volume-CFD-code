!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: The Main Transition subroutine for 3D flows                             //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar, M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir            //!
!// Doc ID: MC5F056F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KWSST_Transition_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,&
                                    NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,&
									DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut,Kinf,Oinf,Ginf)
 Implicit None
!******************************************************************************************* 
 Intent(In   )::Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Mu0,WB,WNP1,Mu,Kinf,&
				Oinf,Ginf
 Intent(InOut  )::WTNP1,Mut

 Integer::Dim,I,J,Ii,Jj,NC,NF,NFW1,NFW2,NF1,NF2,NS,Ncyc,NP,NRKS,ME,NFI1,NFI2,NFO1,NFO2,NFS1,&
          NFS2,NFF1,NFF2
 Real(8)::Mu0,RKco,K,Omega,MR,Co,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,Rho,&
          Kinf,Oinf,Ginf
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:3,1:Dim)::WTNP1,WTB,WTN,Cont,Dift,St
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,DKX_C,DKY_C,&
                           DKZ_C,&
                           DOmegX_C,DOmegY_C,DOmegZ_C,&
                           DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,DOmegZ_F,DGamaX_F,DGamaY_F,&
						   DGamaZ_F,F11,F22,Sigk,Sigw,Sigg,Beta,Gama,X,Y,Z,NX,NY,NZ,Vol,Mu,&
						   Mut,Mutq,XC,YC,ZC,DW,DT,Error
 Real(8),Dimension(1:Dim)::Rev,Strain,Vor,Fonset,Fturb
 Real(8),Dimension(1:5)::RKJ
 Integer,Dimension(1:Dim)::FaceType
 Real(8)::Flength,ce2,ca2,sigmaG,ck,Csep,Relimtc

!*******************************************************************************************
!Part 1:
 Bstar=0.09
 Sigk1=0.85   ;   Sigw1=0.5     ;   Beta1=0.075    ;   Gama1=5.0/9.0
 Sigk2=1.0    ;   Sigw2=0.856   ;   Beta2=0.0828   ;   Gama2=0.44

 Flength=100.0
 ce2 = 50.0
 ca2 = 0.06
 sigmaG=1.0
 ck=1.0
 Csep=1.0
 Relimtc=1100.0
 
!Part 2:
 Do I=1,NC
    WTN(1,I) = WTNP1(1,I)
    WTN(2,I) = WTNP1(2,I)
    WTN(3,I) = WTNP1(3,I)
    Mutq(I)  = Mut(I)
 End Do

!Part 3:
 Do NS=1,NRKS
          
   !Part 4:  
    RKco=1.0/(NRKS-NS+1.0)
          
   !Part 5:
    Call KwSST_Trans_BC3D(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,&
                          MR,NX,NY,NZ,DW,Mu,WB,WNP1,WTNP1,Kinf,Oinf,Ginf,WTB)

   !Part 6: 
    Call Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,&
	                         DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
 
   !Part 7:
    Call KwSST_Trans_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WTNP1,WTB,WB,&
                                   DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C)
 
   !Part 8:
    Call KwSST_Trans_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,&
	                            Z,XC,YC,ZC,NX,NY,NZ,WNP1,WTNP1,WB,WTB,DKX_F,DKY_F,DKZ_F,&
								DOmegX_F,DOmegY_F,DOmegZ_F,DGamaX_F,DGamaY_F,DGamaZ_F)
 
   !part 9:
    Call KwSST_Trans_Func3D(Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,&
                               DOmegZ_C,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,&
							   F11,F22,Sigk,Sigw,Sigg,Beta,Gama,Rev,Strain,Vor,Fonset,Fturb,&
							   DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
   !Part 10:
	Call KwSST_Trans_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)
       
   !Part 11:
	Call KwSST_Trans_Dif3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,&
                            DOmegX_F,DOmegY_F,DOmegZ_F,DGamaX_F,DGamaY_F,DGamaZ_F,&
                            MR,Sigk,Sigw,Sigg,WNP1,WTNP1,Mu,Mut,Dift)       
   !Part 12:
	Call KwSST_Trans_Source3D(Dim,NC,MR,Vol,WNP1,WTNP1,Mut,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,&
                                 DVZ_C,DWX_C,DWY_C,DWZ_C,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,&
								 DomegZ_C,Flength,ce2,ca2,sigmaG,ck,Csep,Relimtc,Bstar,Sigw2,&
								 Beta,Gama,F11,Rev,Strain,Vor,Fonset,Fturb,St,Mu)
 
   !Part 13:
    Do J=1,NC
 
	   Co = RKco*DT(J)/Vol(J)
 
       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) + St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) + St(2,J) ) 
       WTNP1(3,J) = WTN(3,J) - Co*( Cont(3,J) + Dift(3,J) + St(3,J) ) 

      !Part 14: 
       if(WTNP1(1,J)<0.0 )                     WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 )                     WTNP1(2,J)=WTN(2,J)
       if(WTNP1(3,J)<0.0 .OR. WTNP1(3,J)>1.0 ) WTNP1(3,J)=WTN(3,J)

      !Part 15:
       Rho   = WNP1(1,J)
       K     = WTNP1(1,J)/Rho
       Omega = WTNP1(2,J)/Rho

      !Part 16: 
	   Mut(J) = 0.31*Rho*K / max(0.31*Omega,Strain(J)*F22(J)) /MR
	  
    End Do 
   
 END DO

!Part 17:
 Do I=1,NC
	Error(I) = Dabs(Mut(I)-Mutq(I))
 End Do
      
 if ( Mod(Ncyc,10)==0 ) then
  print*, "Error=" , maxval(Error(1:NC))
 end if
!*******************************************************************************************
 End
!###########################################################################################
