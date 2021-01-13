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
 Subroutine KeChien_Source(Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY,St)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY
 Intent(Out  )::St

 Integer::Dim,I,II,NC,ME,P1,P2
 Real(8)::K,Epsilon,Rho,MR,Yn,Ce1,Ce2,Pk,Pe,Txx,Txy,Tyy,Lk,Le,fe1,fe2,Tauwall,Ustar,Yplus,Rt
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
!********************************************************************************************* 
!Part 1:
 Do I=1,NC
           
   !Part 2:
    Rho     = WNP1(1,I)
    k       = WTNP1(1,I)/Rho
    Epsilon = WTNP1(2,I)/Rho
           
   !Part 3:
    II = INW(I)     ! II: Wall Face
    Yn = DW(I) 
    ME = IDS(1,II)

   !Part 4:
    Tauwall = Mu(ME)*DUY(II)
    Ustar   = Dsqrt(abs(MR*tauwall/Wnp1(1,ME)))
    Yplus   = (1.0/MR)*Rho*ustar*Yn/Mu(I)
           
   !Part 5: 
    Rt = (Rho*k*k)/(Mu(I)*Epsilon) / MR
    fe1  =  1.0
    fe2  =  1.0 - (0.4/1.8)*exp(-Rt*Rt/36.0)
 
   !Part 6:
    Txx = MR * Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) - Rho*K/1.5
    Txy = MR * Mut(I)*( DDUY(I)+DDVX(I)  )
    Tyy = MR * Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) - Rho*K/1.5

    Pe =  txx*DDUX(I) + txy*(DDUY(I)+DDVX(I)) + tyy*DDVY(I) 
    Pk = min (Pe,10.0*Rho*Epsilon)
 
   !Part 7:  
    Lk = MR * (2.0*Mu(I)*k)/(Yn*Yn)
    Le = MR * (exp(-0.5*Yplus))*(2.0*Mu(I)*Epsilon)/(Yn*Yn)

    St(1,I) = A(I) * ( -Pk + Rho*Epsilon + Lk ) 
    St(2,I) = A(I) * ( -Ce1*fe1*Epsilon*Pk/k + Ce2*fe2*Rho*Epsilon*Epsilon/k + Le )   

 END DO
!*********************************************************************************************
 End
!###########################################################################################

