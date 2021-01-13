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
 Subroutine KwWilcox_Source3D(Dim,NC,MR,ALFA,BETA,BETA_S,Mut,Vol,WTNP1,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR,ALFA,BETA,BETA_S,Mut,Vol,WTNP1,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Intent(Out  )::St

 Integer::I
 Real(8)::Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
 Real(8)::K,Mum,Uii,Rk,Coeff,Omega,ALFA,BETA,BETA_S,Pk,Rho
!-------------------------------------------------------------------------------------------
 Integer::Dim
 Integer::NC
 Real(8)::MR
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::St
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:2,1:Dim)::WTB
 Real(8),Dimension(1:Dim)::Mut,Vol,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
!*********************************************************************************************
!Part 1:
 Do I=1,NC
      
   !Part 2: 
    Rho   = WNP1(1,I)
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho
    
   !Part 3:
    Uii = (DUX_C(I)+DVY_C(I)+DWZ_C(I))/1.5

    Sxx =2*DUX_C(I)-Uii  ;  Sxy =  DUY_C(I)+DVX_C(I)  ;  Sxz =  DUZ_C(I)+DWX_C(I)
                            Syy =2*DVY_C(I)-Uii       ;  Syz =  DVZ_C(I)+DWY_C(I)
                                                         Szz =2*DWZ_C(I)-Uii
    !Part 4:
	Coeff = MR * Mut(I)
    Rk    = -(2.0/3.0)*Rho*K
    
    Txx = Coeff * Sxx  + Rk  ;  Txy = Coeff * Sxy        ;  Txz = Coeff * Sxz
    Tyx = Txy                ;  Tyy = Coeff * Syy  + Rk  ;  Tyz = Coeff * Syz
    Tzx = Txz                ;  Tzy = Tyz                ;  Tzz = Coeff * Szz   + Rk

    Pk =  Txx*DUX_C(I) + Tyy*DVY_C(I) + Tzz*DWZ_C(I) + Txy*(DUY_C(I)+DVX_C(I)) + Txz*(DWX_C(I)+DUZ_C(I)) + Tyz*(DVZ_C(I)+DWY_C(I))

   !Part 5:
    Pk       =  min (Pk,20.0*BETA_S*Rho*Omega*K)
 
   !Part 6: 
    St(1,I) = Vol(I) * ( Pk - Rho*BETA_S*K*Omega ) 
    St(2,I) = Vol(I) * ( ALFA * PK * Omega/(K+1.e-20) -BETA* Rho * Omega*Omega) 
		
 END DO
!*********************************************************************************************
 End
!###########################################################################################


