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
 Subroutine KwSST_Func3D(Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C,Sigk1,Sigk2,&
                       Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,F11,F22,Sigk,Sigw,Beta,Gama)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C,Sigk1,Sigk2,Sigw1,&
                Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar
 Intent(Out  )::F11,F22,Sigk,Sigw,Beta,Gama

 Integer::Dim,I,NC
 Real(8)::CDKw,F1,F2,Arg1,Arg2,Part1,Part2,k,Omega,Rho,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,&
          Beta2,Gama1,Gama2,Bstar,MR
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::DW,Mu,DKX_C,DKY_C,DKZ_C,DOmegX_C,DOmegY_C,DOmegZ_C,F11,F22,Sigk,Sigw,Beta,Gama
!*********************************************************************************************

!Part 1:
 Do I=1,NC

   !Part 2:
    Rho   = WNP1(1,I)
    k     = WTNP1(1,I)/Rho
    Omega = WTNP1(2,I)/Rho

   !Part 3:    
    CDkw = max( 2.0*Rho*Sigw2*(DKX_C(I)*DomegX_C(I)+DKY_C(I)*DomegY_C(I)+DKZ_C(I)*DomegZ_C(I))/Omega,1.0D-20)

   !Part 4:
	Part1 = Dsqrt(K) / ( Bstar*Omega*DW(I) )
	Part2 = 500.0*Mu(I)*MR / (DW(I)*DW(I)*Rho*Omega)     

    arg1=min ( max(Part1,Part2), 4.0*Rho*Sigw2*K / (DW(I)*DW(I)*CDkw) )
    arg2=max(2*Part1,Part2)

    F1=tanh(arg1*arg1*arg1*arg1)
    F2=tanh(arg2*arg2)

   !Part 5:
	F11(I) = F1
	F22(I) = F2

   !Part 6:
    Sigk(I) = F1 * Sigk1 + (1.0-F1) * Sigk2
    Sigw(I) = F1 * Sigw1 + (1.0-F1) * Sigw2
    Beta(I) = F1 * Beta1 + (1.0-F1) * Beta2
    Gama(I) = F1 * Gama1 + (1.0-F1) * Gama2
  
 End Do
!*********************************************************************************************
 End
!###########################################################################################

