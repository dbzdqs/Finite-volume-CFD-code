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
 Subroutine KwSST_Source3D(Dim,NC,MR,Vol,WNP1,WTNP1,Mut,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,&
                         DWX_C,DWY_C,DWZ_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C,Bstar,Sigw2,Beta,Gama,F11,St)
 
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR,Vol,WNP1,WTNP1,Mut,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,&
                DWX_C,DWY_C,DWZ_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C,Bstar,Sigw2,Beta,Gama,F11
 Intent(Out  )::St

 Integer::Dim,I,NC
 Real(8)::MR,Bstar,Sigw2,K,Omega,Rho,Sxx,Syy,Szz,Sxy,Sxz,Syz,Rk,Coeff,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz,Pk,Uii
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::Mut,Beta,Gama,F11,Vol,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,DWY_C,DWZ_C ,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
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

	Coeff = MR * Mut(I)
    Rk    = -(2.0/3.0)*Rho*K
    
    Txx = Coeff * Sxx  + Rk  ;  Txy = Coeff * Sxy        ;  Txz = Coeff * Sxz
    Tyx = Txy                ;  Tyy = Coeff * Syy  + Rk  ;  Tyz = Coeff * Syz
    Tzx = Txz                ;  Tzy = Tyz                ;  Tzz = Coeff * Szz   + Rk
    
 !   Txx = MR * Mut(I) * ( 2*DUX_C(I) - Uii      ) -(2.0/3.0)*Rho*K
	!Txy = MR * Mut(I) * (   DVX_C(I) + DUY_C(I) )
 !   Tyy = MR * Mut(I) * ( 2*DVY_C(I) - Uii      ) -(2.0/3.0)*Rho*K

   !Part 4:
    Pk = Txx*DUX_C(I) + Tyy*DVY_C(I) + Tzz*DWZ_C(I) + Txy*(DUY_C(I)+DVX_C(I)) + Txz*(DWX_C(I)+DUZ_C(I)) + Tyz*(DVZ_C(I)+DWY_C(I))

   !Part 5:
    Pk = min (Pk,20.0*Bstar*Rho*Omega*K)

   !Part 6:
    St(1,I) = Vol(I) * ( Pk - Bstar*Rho*Omega*K  ) 
    St(2,I) = Vol(I) * ( (Gama(I)*Rho/Mut(I)/MR)*Pk - Beta(I)*Rho*Omega*Omega &
              + 2.0*(1.0-F11(I))*Rho*Sigw2*(DKX_C(I)*DomegX_C(I)+DKY_C(I)*DomegY_C(I)+DKZ_C(I)*DomegZ_C(I))/Omega  )    
	  	   		  
 END DO
!********************************************************************************************* 
 End
!###########################################################################################

