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
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KwSST_Sust_Source(Dim,NC,MR,A,WNP1,WTNP1,Mut,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,&
	                     DOmegX_C,DOmegY_C,Bstar,Sigw2,Beta,Gama,F11,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR,A,WNP1,WTNP1,Mut,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DOmegX_C,&
	            DOmegY_C,Bstar,Sigw2,Beta,Gama,F11
 Intent(Out  )::St

 Integer::Dim,I,NC
 Real(8)::MR,Bstar,Sigw2,K,Omega,Rho,Txx,Txy,Tyy,Pk,Uii
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::Mut,Beta,Gama,F11,A,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
!********************************************************************************************* 
!Part 1:
 Do I=1,NC
     
   !Part 2: 
    Rho   = WNP1(1,I)
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho
    
   !Part 3:
    Uii = (DUX_C(I)+DVY_C(I))/1.5

    Txx = MR * Mut(I) * ( 2*DUX_C(I) - Uii      ) -(2.0/3.0)*Rho*K
	Txy = MR * Mut(I) * (   DVX_C(I) + DUY_C(I) )
    Tyy = MR * Mut(I) * ( 2*DVY_C(I) - Uii      ) -(2.0/3.0)*Rho*K

   !Part 4:
    Pk = Txx*DUX_C(I) + Txy*(DUY_C(I)+DVX_C(I)) + Tyy*DVY_C(I)

   !Part 5:
    Pk = min (Pk,20.0*Bstar*Rho*Omega*K)

   !Part 6:
    St(1,I) = A(I) * ( Pk - Bstar*Rho*Omega*K  ) 
    St(2,I) = A(I) * ( (Gama(I)*Rho/Mut(I)/MR)*Pk - Beta(I)*Rho*Omega*Omega &
              + 2.0*(1.0-F11(I))*Rho*Sigw2*(DKX_C(I)*DomegX_C(I)+DKY_C(I)*DomegY_C(I))/Omega  )    
	  	   		  
 END DO
!********************************************************************************************* 
 End
!###########################################################################################

