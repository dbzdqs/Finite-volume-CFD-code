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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine SA_ProdDest3D(Dim,NC,Vol,Dw,MR,Cv1,Cb1,Cw1,Cw2,Cw3,Kei,WNP1,Mu,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,DWY_C,DWZ_C,WTNP1,Dest,Prod)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,Vol,Dw,MR,Cv1,Cb1,Cw1,Cw2,Cw3,Kei,WNP1,Mu,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,DWY_C,DWZ_C,WTNP1
 Intent(Out  )::Dest,Prod

 Integer::Dim,I,NC
 Real(8)::Shat,Nu,Fv1,Fv2,Chi,Chi3,S,K2D2,Kei,R,G,Cv13,Cw36,Cv1,Cw1,Cw2,Cw3,T2,Fw,MR,Cb1,D2,Vor,Vorticity,StrainRateMag,Ft2
 Real(8),Dimension(1:Dim)::Vol,Dw,Mu
 Real(8)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C ,DWX_C,DWY_C,DWZ_C 
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:2,1:Dim)::WTNP1,Dest,Prod
!*********************************************************************************************
!Part 1:
 Cv13=Cv1*Cv1*Cv1
 Cw36=Cw3**6.0

!Part 2:
 Do I=1,NC

   !Part 3:
	Chi  = WTNP1(1,I)/Mu(I)
	Chi3 = Chi*Chi*Chi

   !Part 4:
	Fv1 = Chi3/(Chi3+Cv13)

   !Part 5:
	Fv2 = 1.0 - Chi/(1+Chi*Fv1)
    
    DUX = DUX_C(I) ; DUY = DUY_C(I) ; DUZ = DUZ_C(I)
    DVX = DVX_C(I) ; DVY = DVY_C(I) ; DVZ = DVZ_C(I)
    DWX = DWX_C(I) ; DWY = DWY_C(I) ; DWZ = DWZ_C(I)

   !Part 6:
	Vor = Vorticity(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ)
	S   = StrainRateMag(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ) !+  2*Min(0.0 , Vor-Omega)

   !Part 7:
	D2   = Dw(I)*Dw(I)
	K2D2 = Kei*Kei*D2

   !Part 8:
    Shat = max(S + MR * ( WTNP1(1,I)/(K2D2*WNP1(1,I)) )*Fv2,0.0)

   !Part 9:
    Ft2=1.2*exp(-0.5*Chi*Chi)
    Prod(1,I) = Vol(I) * Cb1 * (1-Ft2) * WTNP1(1,I) * Shat

   !Part 10:
    R  = MR *  WTNP1(1,I) / (WNP1(1,I)*Shat*K2D2)
	if(R>10.) R=10.

	G  = R+Cw2*(R**6. - R)
	Fw = G*( (1+Cw36) / (G**6. + Cw36) )**(1/6.)

   !Part 11:
    Dest(1,I) = MR * Vol(I) * (Cw1*Fw-(Cb1/Kei*Kei)*Ft2)  *  WTNP1(1,I)*WTNP1(1,I) / (WNP1(1,I)*D2)
 End do
!*********************************************************************************************
 End
!###########################################################################################

