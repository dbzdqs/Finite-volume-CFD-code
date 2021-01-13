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
 Subroutine KwBredberg_Source(Dim,NC,MR,CK,CMU,CW,CW1,CW2,Mut,A,WTNP1,WNP1,DUX_C,DUY_C,DVX_C,&
                              DVY_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR,CK,CMU,CW,CW1,CW2,Mut,A,WTNP1,WNP1,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,&
                DKY_C,DOmegX_C,DOmegY_C
 Intent(Out  )::St

 Integer::Dim,I,ME,NE,NC,NF,NF1,NF2
 Real(8)::K,Mum,Txx,Tyy,Txy,MR
 Real(8)::Omega,Uii,RHS3
 Real(8)::CK,CMU,CW,CW1,CW2,Pk,Rho
 Real(8)::CR1,CR2,CR3,W11,W12,W22,S11,S12,S22,Cscale,NORM_SW
 Real(8)::Rhat,Rstar,Frot,Fr1,Fr,WW,SS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:2,1:Dim)::St,WTNP1,WTB
 Real(8),Dimension(1:Dim)::Mut,DUX_C,DUY_C,DVX_C,DVY_C
 Real(8),Dimension(1:Dim)::DKX_C,DKY_C,DOmegX_C,DOmegY_C
 Real(8),Dimension(1:Dim)::A
!*********************************************************************************************
!Part 1:
 CR1=1.0
 CR2=2.0
 CR3=1.0
 Cscale=3.0

!Part 2:
 Do I=1,NC
      
   !Part 3: 
    Rho   = WNP1(1,I)
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho
    
   !Part 4:
    Uii = (DUX_C(I)+DVY_C(I))/1.5

    Txx = MR * Mut(I) * ( 2*DUX_C(I) - Uii      ) -(2.0/3.0)*Rho*K
	Txy = MR * Mut(I) * (   DVX_C(I) + DUY_C(I) )
    Tyy = MR * Mut(I) * ( 2*DVY_C(I) - Uii      ) -(2.0/3.0)*Rho*K

   !Part 5:
    Pk = Txx*DUX_C(I) + Txy*(DUY_C(I)+DVX_C(I)) + Tyy*DVY_C(I)

   !Part 6:
	S11  = DUX_C(I)
	S12  = 0.50*(DUY_C(I)+DVX_C(I)) 
    S22  = DVY_C(I)
    
	!Part 7: 
    W11  =  1.e-30
	W12  =  0.50*(DUY_C(I)-DVX_C(I))
	W22  =  1.e-30 
	
   !Part 8:
	WW   = W11*W11+2.0*W12*W12+W22*W22
	SS   = S11*S11+2.0*S12*S12+S22*S22
	NORM_SW=(SQRT(2.0*SS))/(SQRT(2.0*WW))

    Rhat =(1.0/NORM_SW ) *(1.0/NORM_SW-1.0)  
	Rstar  = NORM_SW
	
   !Part 9:
    Frot   = (1.0+ CR1) * ((2.0*Rstar)/(1.0+Rstar) ) * &
	         (1.0-CR3*ATAN(CR2*Rhat)) - CR1
    !Part 10:
    Fr1    = MAX(MIN(Frot,1.250) , 0.0)
    Fr     = Max(0.0 ,1.0+ Cscale*(Fr1-1.0))
 
   !Part 11:
    Pk       =  min (Fr*Pk,20.0*CK*Rho*Omega*K)
    
   !Part 12:
    RHS3     = (CW/(K+1.e-20))*(Mut(I))*(DKX_C(I)*DOmegX_C(I) + DKY_C(I)*DOmegY_C(I))*MR
   
   !Part 13: 
    St(1,I) = A(I) * ( Pk - Rho*CK*K*Omega ) 
    St(2,I) = A(I) * ( CW1 * PK * Omega/(K+1.e-20) -CW2* Rho * Omega*Omega + Rho*RHS3) 
    
 END DO
!*********************************************************************************************
 End
!###########################################################################################


