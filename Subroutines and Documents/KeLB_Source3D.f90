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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeLB_Source3D(Dim,NC,IDS,DW,Vol,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,GM,P,WTNP1,&
                        Mu,Mut,DUY,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,IDS,DW,Vol,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,GM,P,WTNP1,&
                Mu,Mut,DUY,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Intent(Out  )::St

 Integer::Dim,I,II,NC,ME,P1,P2
 Real(8)::K,MR,Tauwall,Ustar,Yplus
 Real(8)::Lam1,Lam2,Lam3,PKS,CC,GM,Mtu2,Epsilon,PD,TComp,EPS0,F11,Ret,Rets
 Real(8)::Ceps1,Ceps2,Pe,Pk,Yn,Rho,F22,FMU
 Real(8)::CR1,CR2,CR3,W11,W12,W22,S11,S12,S22,Cscale,NORM_SW
 Real(8)::Rhat,Rstar,Frot,Fr1,Fr
 Real(8)::Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz,Wxx,Wxy,Wxz,Wyx,Wyy,Wyz,Wzx,Wzy,Wzz,&
           SijSij,WijWij,Uii,Rk,Coeff,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::Vol
 Real(8),Dimension(1:Dim)::DW
 Real(8),Dimension(1:Dim)::Mu
 Real(8),Dimension(1:Dim)::Mut
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Real(8),Dimension(1:Dim)::DUY
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:2,1:Dim)::St
!*********************************************************************************************
!part 1:
 CR1=1.0
 CR2=2.0
 CR3=1.0
 Cscale=1.0
  
!Part 2:
 Do I=1,NC
           
   !Part 3:
    Rho     = WNP1(1,I)
    k       = WTNP1(1,I)/Rho
    Epsilon = WTNP1(2,I)/Rho
           
   !Part 4:
    II = INW(I)     ! II: Wall Face
    Yn = DW(I) 
    ME = IDS(1,II)

   !Part 5:
    Tauwall = Mu(ME)*DUY(II) !!!!!!!!!!!!!!!!!!!!!!!!!
    Ustar   = Dsqrt(abs(MR*tauwall/Wnp1(1,ME)))
    Yplus   = (1.0/MR)*Rho*ustar*Yn/Mu(I)
           
   !Part 6: 
    Ret  = (K * K * Rho)/(Mu(I)*Epsilon)
    Rets = (1.0/MR)* Ret
    
   !Part 7:
    FMu = 0.04 + (1.0 - 0.04)*((1.0 - EXP(-(YPLUS - 8.0)/26.0))**2)
    F11 = 1.0   + (0.05/(FMu+1.E-10))**3.0  	  
    F22 = (1.0  -0.22* EXP(-Rets * Rets-1.E-10) ) 
       
   !Part 8:
    Uii = (DUX_C(I)+DVY_C(I)+DWZ_C(I))/1.5

    Sxx =2*DUX_C(I)-Uii  ;  Sxy =  DUY_C(I)+DVX_C(I)  ;  Sxz =  DUZ_C(I)+DWX_C(I)
                            Syy =2*DVY_C(I)-Uii       ;  Syz =  DVZ_C(I)+DWY_C(I)
                                                         Szz =2*DWZ_C(I)-Uii

	Coeff = MR * Mut(I)
    Rk    = -(2.0/3.0)*Rho*K
    
    Txx = Coeff * Sxx  + Rk  ;  Txy = Coeff * Sxy        ;  Txz = Coeff * Sxz
    Tyx = Txy                ;  Tyy = Coeff * Syy  + Rk  ;  Tyz = Coeff * Syz
    Tzx = Txz                ;  Tzy = Tyz                ;  Tzz = Coeff * Szz   + Rk

    Pe =  Txx*DUX_C(I) + Tyy*DVY_C(I) + Tzz*DWZ_C(I) + Txy*(DUY_C(I)+DVX_C(I)) + Txz*(DWX_C(I)+DUZ_C(I)) + Tyz*(DVZ_C(I)+DWY_C(I))

   !Part 9:
    Sxx = 0.5 * (DUX_C(I)+DUX_C(I)) 
    Sxy = 0.5 * (DUY_C(I)+DVX_C(I))
    Sxz = 0.5 * (DUZ_C(I)+DWX_C(I))
    
    Syx = 0.5 * (DVX_C(I)+DUY_C(I))
    Syy = 0.5 * (DVY_C(I)+DVY_C(I))
    Syz = 0.5 * (DVZ_C(I)+DWY_C(I))
    
    Szx = 0.5 * (DWX_C(I)+DUZ_C(I))
    Szy = 0.5 * (DWY_C(I)+DVZ_C(I))
    Szz = 0.5 * (DWZ_C(I)+DWZ_C(I))
   
    !Part 10:
    Wxx = 0.5 * (DUX_C(I)-DUX_C(I))
    Wxy = 0.5 * (DUY_C(I)-DVX_C(I))
    Wxz = 0.5 * (DUZ_C(I)-DWX_C(I))
    
    Wyx = 0.5 * (DVX_C(I)-DUY_C(I))
    Wyy = 0.5 * (DVY_C(I)-DVY_C(I))
    Wyz = 0.5 * (DVZ_C(I)-DWY_C(I))
    
    Wzx = 0.5 * (DWX_C(I)-DUZ_C(I))
    Wzy = 0.5 * (DWY_C(I)-DVZ_C(I))
    Wzz = 0.5 * (DWZ_C(I)-DWZ_C(I))
    
    !Part 11:
    SijSij = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + &
             Syx*Syx + Syy*Syy + Syz*Syz + &
             Szx*Szx + Szy*Szy + Szz*Szz

    WijWij = Wxx*Wxx + Wxy*Wxy + Wxz*Wxz + &
             Wyx*Wyx + Wyy*Wyy + Wyz*Wyz + &
             Wzx*Wzx + Wzy*Wzy + Wzz*Wzz
      
	NORM_SW=(SQRT(2.0*SijSij))/(SQRT(2.0*WijWij))

    Rhat =(1.0/NORM_SW ) *(1.0/NORM_SW-1.0)  	
	Rstar  = NORM_SW

   !Part 12:
    Frot   = (1.0+ CR1) * ((2.0*Rstar)/(1.0+Rstar) ) * &
	         (1.0-CR3*ATAN(CR2*Rhat)) - CR1
   !Part 13:
    Fr1    = MAX(MIN(Frot,1.25d0) , 0.0)
    Fr     = Max(0.0 ,1.0+ Cscale*(Fr1-1.0))

   !Part 14:
    Pk       =   min (Fr*Pe,10.0*Rho*Epsilon)
 
   !Part 15:  
	CC       =  SQRT(GM*P(I)/ Rho)
    Mtu2     =  (2.0*K)/(CC*CC)
    EPS0     =  Lam1 * Mtu2 *  Epsilon
    TComp    =  Rho*(Epsilon+EPS0)
    PD       =  (- Lam2 * PK * Mtu2)+ (Lam3 * Rho * Epsilon * Mtu2  )
 
   !Part 16: 
    St(1,I) = Vol(I) * (- Pk + TComp - PD ) 
    St(2,I) = Vol(I) * ( -Ceps1* F11 * PK * Epsilon/(K+1.e-20) +Ceps2 * F22 * Rho * Epsilon*Epsilon/(K+1.e-20)) 

 END DO
!*********************************************************************************************
 End
!###########################################################################################

