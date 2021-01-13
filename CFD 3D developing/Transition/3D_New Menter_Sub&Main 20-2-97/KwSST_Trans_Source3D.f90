!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate Source Term of Kw New_Menter Transition Model              //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar, M.A.Zolajani Iran, Tehran, OpenFlows@chmail.ir              //!
!// Doc ID: MC2F061F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_Source3D(Dim,NC,MR,Vol,WNP1,WTNP1,Mut,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,&
                                 DVZ_C,DWX_C,DWY_C,DWZ_C,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,&
								 DomegZ_C,Flength,ce2,ca2,sigmaG,ck,Csep,Relimtc,Bstar,Sigw2,&
								 Beta,Gama,F11,Rev,Strain,Vor,Fonset,Fturb,St,Mu)
 
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,MR,Vol,WNP1,WTNP1,Mut,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,&
                DWZ_C,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DomegZ_C,Bstar,Sigw2,Beta,Gama,F11,&
				Rev,Strain,Vor,Fonset,Fturb,Mu
               				
 Intent(Out  )::St

 Integer::Dim,I,NC
 Real(8)::MR,Bstar,Sigw2,K,Omega,Rho,Sxx,Syy,Szz,Sxy,Sxz,Syz,Rk,Coeff,Txx,Txy,Txz,Tyx,Tyy,Tyz,&
          Tzx,Tzy,Tzz,Pk,Uii,DK,Pbark
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::Mut,Beta,Gama,F11,Vol,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,&
                           DWY_C,DWZ_C,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C
 Real(8),Dimension(1:3,1:Dim)::WTNP1,St
 Real(8)::Flength,ce2,ca2,sigmaG,ck,Csep,Relimtc,Gamma,DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ
 Real(8)::Pw,Modi,Flimon,Pgama,Egama,Plimk,Dbark,StrainRateMag,Vorticity
 Real(8),Dimension(1:Dim)::Rev,Mu,Vor,Fturb,Strain,Fonset
!******************************************************************************************* 
!Part 1:
 Flength=100.0
 ce2 = 50.0
 ca2 = 0.06
 sigmaG=1.0
 ck=1.0
 Csep=1.0
 Relimtc=1100.0

!Part 2:
Do I=1,NC
     St(1,I)=0.0
     St(2,I)=0.0
     St(3,I)=0.0
 End Do

 !Part 3:
 Do I=1,NC
     
    Rho   = WNP1(1,I)
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho
    Gamma = WTNP1(3,I)/Rho
    
    DUX = DUX_C(I) ; DUY = DUY_C(I) ; DUZ = DUZ_C(I)
    DVX = DVX_C(I) ; DVY = DVY_C(I) ; DVZ = DVZ_C(I)
    DWX = DWX_C(I) ; DWY = DWY_C(I) ; DWZ = DWZ_C(I)
    
 

   !Part 4:
	Pk = Mut(I)*Strain(I)*Vor(I)*MR    !kato Launder formula
    
   !Part 5:
    Modi = max(Gamma,0.1)

    Flimon=min(max(Rev(I)/(2.2*Relimtc)-1.0,0.0),3.0)
    Plimk=5.0*ck*max(Gamma-0.2,0.0)*(1.0-Gamma)*Flimon*max(3.0*Csep*Mu(I)*MR-Mut(I)*MR,0.0)*&
	      Strain(I)*Vor(I)
    
    Pgama=Flength*Rho*Strain(I)*Gamma*(1.0-Gamma)*Fonset(I)
    Egama=ca2*Rho*Vor(I)*Gamma*Fturb(I)*(ce2*Gamma-1.0)
    
  
   !Part 6:???????????????????????
    St(1,I) = Vol(I) * ( - Gamma*Pk + Bstar*Rho*Omega*K*Modi - Plimk )
   
    St(2,I) = Vol(I) * ( - (Gama(I)*Rho/Mut(I)/MR)*Pw + Beta(I)*Rho*Omega*Omega - 2.0*(1.0-F11(I))*&
	          Rho*Sigw2*(DKX_C(I)*DomegX_C(I)+DKY_C(I)*DomegY_C(I)+DKZ_C(I)*DomegZ_C(I))/Omega  )
	                        
    St(3,I) = Vol(I) * (- Pgama + Egama)
   

 END DO

!******************************************************************************************* 
 End
!###########################################################################################
