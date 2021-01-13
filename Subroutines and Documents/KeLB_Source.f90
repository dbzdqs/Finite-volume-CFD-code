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
!// Developed by: H. Nazari, Mechanical Eng., Amirkabir University of Technology           //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeLB_Source(Dim,NC,IDS,DW,A,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,GM,P,WTNP1,&
                        Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,IDS,DW,A,INW,MR,Ceps1,Ceps2,WNP1,WTNP1,Mu,Mut,DUY,DDUX,DDUY,DDVX,&
                DDVY,Lam1,Lam2,Lam3,P,GM
 Intent(Out  )::St

 Integer::Dim,I,II,NC,ME,P1,P2
 Real(8)::K,MR,Txx,Txy,Tyy,Tauwall,Ustar,Yplus
 Real(8)::Lam1,Lam2,Lam3,PKS,CC,GM,Mtu2,Epsilon,PD,TComp,EPS0,F11,Ret,Rets
 Real(8)::Ceps1,Ceps2,Pe,Pk,Yn,Rho,F22,FMU
 Real(8)::CR1,CR2,CR3,W11,W12,W22,S11,S12,S22,Cscale,NORM_SW
 Real(8)::Rhat,Rstar,Frot,Fr1,Fr,WW,SS
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY,P
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
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
    Tauwall = Mu(ME)*DUY(II)
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
    Txx = MR * Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) - Rho*K/1.5
    Txy = MR * Mut(I)*( DDUY(I)+DDVX(I)  )
    Tyy = MR * Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) - Rho*K/1.5

    Pe =  txx*DDUX(I) + txy*(DDUY(I)+DDVX(I)) + tyy*DDVY(I) 

   !Part 9:
	S11  = DDUX(I)
	S12  = 0.50*(DDUY(I)+DDVX(I)) 
    S22  = DDVY(I)

   !Part 10:	 
    W11  =  1.e-30
	W12  =  0.50*(DDUY(I)-DDVX(I))
	W22  =  1.e-30

   !Part 11:
	WW   = W11*W11+2.0*W12*W12+W22*W22
	SS   = S11*S11+2.0*S12*S12+S22*S22
	NORM_SW=(SQRT(2.0*SS))/(SQRT(2.0*WW))

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
    St(1,I) = A(I) * (- Pk + TComp - PD ) 
    St(2,I) = A(I) * ( -Ceps1* F11 * PK * Epsilon/(K+1.e-20) +Ceps2 * F22 * Rho * Epsilon*Epsilon/(K+1.e-20)) 

 END DO
!*********************************************************************************************
 End
!###########################################################################################

