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
!// Date: Mar., 05, 2013                                                                   //!
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
 Subroutine KwSST_V_Source(Dim,NC,MR,A,Wnp1,Wntp1,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,&
                         DDOmegY,Bstar,Sigw2,Beta,Gama,F11,St)
 Implicit None
!*********************************************************************************************

 Integer::Dim,I,NC
 Real(8)::MR,Bstar,Sigw2,K,Omega,Rho,Txx,Txy,Tyy,Pk,Pw,Uii,Vor
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,Beta,Gama,F11,A
 Real(8),Dimension(1:2,1:Dim)::Wntp1,St
!********************************************************************************************* 
!Part 2:
 Do I=1,NC
     
   !Part 3: 
    Rho   = Wnp1(1,I)
    k     = Wntp1(1,I)/Rho
    Omega = Wntp1(2,I)/Rho
    
   !Part 4:
    Uii = (DDUX(I)+DDVY(I))/1.5

    Txx = MR * Mut(I) * ( 2*DDUX(I) - Uii     ) -(2.0/3.0)*Rho*K
	Txy = MR * Mut(I) * (   DDVX(I) + DDUY(I) )
    Tyy = MR * Mut(I) * ( 2*DDVY(I) - Uii     ) -(2.0/3.0)*Rho*K
    
    Vor =  Dsqrt ( (DDUY(I)-DDVX(I))*(DDUY(I)-DDVX(I)) )

    !Part 5:
     Pw = MR * Mut(I)*Vor*Vor-(2.0/3.0)*Rho*k*(DDUX(I)+DDVY(I)) 
     Pw = abs(Pw)
    
    !Part 6:
    Pk = min (Pw,20.0*Bstar*Rho*Omega*K)

    !Part 7:
    St(1,I) = A(I) * ( Pk - Bstar*Rho*Omega*K  )  
    St(2,I) = A(I) * ( (Gama(I)*Rho/Mut(I)/MR)*Pw - Beta(I)*Rho*Omega*Omega &
              + 2.0*(1.0-F11(I))*Rho*Sigw2*(DDKX(I)*DDomegX(I)+DDKY(I)*DDomegY(I))/Omega  )   
  
 END DO

!********************************************************************************************* 
 End
!###########################################################################################

