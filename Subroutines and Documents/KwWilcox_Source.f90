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
 Subroutine KwWilcox_Source(Dim,NC,MR,ALFA,BETA,BETA_S,Mut,A,WTNP1,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR,Mut,WNP1,WB,ALFA,BETA,BETA_S,WTNP1,DUX_C,DUY_C,DVX_C,DVY_C
 Intent(Out  )::St

 Integer::Dim,I,NC
 Real(8)::K,Mum,Txx,Tyy,Txy,MR,Omega,Uii,ALFA,BETA,BETA_S,Pk,Rho
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::St,WTNP1,WTB
 Real(8),Dimension(1:Dim)::Mut,A,DUX_C,DUY_C,DVX_C,DVY_C
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
    Pk       =  min (Pk,20.0*BETA_S*Rho*Omega*K)
 
   !Part 6: 
    St(1,I) = A(I) * ( Pk - Rho*BETA_S*K*Omega ) 
    St(2,I) = A(I) * ( ALFA * PK * Omega/(K+1.e-20) -BETA* Rho * Omega*Omega) 
		
 END DO
!*********************************************************************************************
 End
!###########################################################################################


