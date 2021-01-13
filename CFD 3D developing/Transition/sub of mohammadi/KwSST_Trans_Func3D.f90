!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate Functions of Transition Model                              //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar, M.A.Zoljanahi Iran, Tehran, OpenFlows@chmail.ir             //!
!// Doc ID: MC2F059F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_Func3D(Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,&
                               DOmegZ_C,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,&
							   F11,F22,Sigk,Sigw,Sigg,Beta,Gama,Rev,Strain,Vor,Fonset,Fturb,&
							   DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C,Sigk1,&
                Sigk2,Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,&
				DVZ_C ,DWX_C,DWY_C,DWZ_C
 Intent(Out  )::F11,F22,Sigk,Sigw,Sigg,Beta,Gama,Rev,Strain,Vor,Fonset,Fturb

 Integer::Dim,I,NC
 Real(8)::CDKw,F1,F2,Arg1,Arg2,Part1,Part2,k,Omega,Rho,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,&
          Beta2,Gama1,Gama2,Bstar,MR,Ry,F3,DUDS,UU
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:3,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::DW,Mu,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C,F11,F22,Sigk,Sigw,&
                           Beta,Gama
 Real(8)::U,V,W,Gamma,DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C 
 Real(8),Dimension(1:Dim)::Fonset,Rev,Strain,Fturb,Vor,Sigg
 Real(8)::CPG1,CPG2,CPG3,ClimPG1,CTU1,CTU2,CTU3,ClimPG2,FPG,Ltl,TUL,Retc,Fonset1,Fonset2,Fonset3,&
          RT,StrainRateMag,Vorticity
!*******************************************************************************************
 
 !Part 1:
 CTU1=100.0
 CTU2=1000.0
 CTU3=1.0
 
 CPG1=14.68    
 ClimPG1=1.5
 CPG2=-7.34
 CPG3=0.0
 ClimPG2=3.0


 Do I=1,NC

   !Part 2:
    Rho   = WNP1(1,I)
    U     = WNP1(2,I)/Rho
    V     = WNP1(3,I)/Rho
    W     = WNP1(4,I)/Rho
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho
    Gamma = WTNP1(3,I)/Rho

    DUX = DUX_C(I) ; DUY = DUY_C(I) ; DUZ = DUZ_C(I)
    DVX = DVX_C(I) ; DVY = DVY_C(I) ; DVZ = DVZ_C(I)
    DWX = DWX_C(I) ; DWY = DWY_C(I) ; DWZ = DWZ_C(I)
    
    Strain(I) =  StrainRateMag(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ)
    Vor(I) = Vorticity(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ)
    
   !Part 3:    
    CDkw = max( 2.0*Rho*Sigw2*(DKX_C(I)*DomegX_C(I)+DKY_C(I)*DomegY_C(I)+DKZ_C(I)*DomegZ_C(I))/Omega,10.0**-10)
    
   !Part 4:
	Part1 = Dsqrt(abs(K)) / ( Bstar*Omega*DW(I) )
	Part2 = 500.0*Mu(I)*MR / (DW(I)*DW(I)*Rho*Omega)     

    arg1=min ( max(Part1,Part2), 4.0*Rho*Sigw2*K / (DW(I)*DW(I)*CDkw) )
    arg2=max(2*Part1,Part2)

    F1=tanh(arg1*arg1*arg1*arg1)
    F2=tanh(arg2*arg2)

    Ry = Rho*DW(I)*dsqrt(abs(K))/(Mu(I)*MR)
    F3 = exp(-(Ry/120.0)**8.0)
    F1 = max( F1,F3 )
    
   !Part 5:
	F11(I) = F1
	F22(I) = F2

   !Part 6:
    Sigk(I) = F1 * Sigk1 + (1.0-F1) * Sigk2
    Sigw(I) = F1 * Sigw1 + (1.0-F1) * Sigw2
    Beta(I) = F1 * Beta1 + (1.0-F1) * Beta2
    Gama(I) = F1 * Gama1 + (1.0-F1) * Gama2
  
    
   !Part 7:
	UU = Dsqrt(U*U + V*V + W*W)
    DUDS =  ((u*u)/(UU*UU)) * DUX +  ((v*v)/(UU*UU)) * DVY + ((W*W)/(UU*UU)) * DWZ + &
	        ((u*v)/(UU*UU)) * (DVX + DUY) + ((u*W)/(UU*UU)) * (DUZ + DWX) + ((V*W)/(UU*UU)) * (DVZ + DWY)

    Ltl=7.57*10**(-3.0)*DUDS*Rho*Dw(I)*Dw(I)/(Mu(I)*MR)+0.0128

    Ltl=min(max(Ltl,-1.0),1.0)
    
    if (Ltl>=0.0) then
     FPG=min(1+CPG1*Ltl,ClimPG1)
    else
     FPG=min(1+CPG2*Ltl+CPG3*min(Ltl+0.0681,0.0),ClimPG2)
    end if
    
    FPG=max(FPG,0.0)
    
    TUL=min(100.0*Dsqrt(2.0*abs(k)/3.0)/(Omega*Dw(I)),100.0)
    Retc=CTU1+CTU2*exp( -CTU3*TUL*FPG )

     
    !Part 8:
    Rev(I) = Rho*Strain(I)*DW(I)*DW(I)/ (Mu(I)*MR)    
    Fonset1 = Rev(I) / (2.2*Retc)   
    Fonset2 = min(Fonset1,2.0)

    RT = Rho*k / (MR*Mu(I)*Omega)
    Fonset3 = max( 1.0 - (RT/3.5)**3.0 , 0.0)
    
    Fonset(I) = max(Fonset2-Fonset3 , 0.0)

    Fturb(I) = exp(-(RT/2.0)**4.0)
    
 End Do
 Sigg(:)=1.0
!*******************************************************************************************
    End
!###########################################################################################
