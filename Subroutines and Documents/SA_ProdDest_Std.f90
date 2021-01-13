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
 Subroutine SA_ProdDest_Std(Dim,NC,A,Dw,MR,Cv1,Cb1,Cw1,Cw2,Cw3,Kei,WNP1,Mu,DUY,DVX,DUX,DVY,WTNP1,Dest,Prod)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,A,Dw,MR,Cv1,Cb1,Cw1,Cw2,Cw3,Kei,WNP1,Mu,DUY,DVX,DUX,DVY,WTNP1
 Intent(Out  )::Dest,Prod

 Integer::Dim,I,NC
 Real(8)::Shat,Nu,Fv1,Fv2,Chi,Chi3,S,K2D2,Kei,R,G,Cv13,Cw36,Cv1,Cw1,Cw2,Cw3,T2,Fw,MR,Cb1,D2,Vor,Omega,Ft2
 Real(8),Dimension(1:Dim)::A,Dw,WTNP1,Mu,DUY,DVX,Dest,Prod,DUX,DVY
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*********************************************************************************************
!Part 1:
 Cv13=Cv1*Cv1*Cv1
 Cw36=Cw3**6.0

!Part 2:
 Do I=1,NC

   !Part 3:
	Chi  = WTNP1(I)/Mu(I)
	Chi3 = Chi*Chi*Chi

   !Part 4:
	Fv1 = Chi3/(Chi3+Cv13)

   !Part 5:
	Fv2 = 1.0 - Chi/(1+Chi*Fv1)

   !Part 6:
    Omega = Dabs( DUY(I)-DVX(I) )
	Vor   = sqrt(2.)*sqrt( DUX(I)*DUX(I) + DVY(I)*DVY(I) + 0.5*(DUY(I)+DVX(I))*(DUY(I)+DVX(I)) )
	S = Omega !+  2*Min(0.0 , Vor-Omega)
   !agar tikeye 2vom S ra dar nazar begirim mishavad vorticity/strain-based vagar na vorticity-based

   !Part 7:
	D2   = Dw(I)*Dw(I)
	K2D2 = Kei*Kei*D2

   !Part 8:
    Shat = max(S + MR * ( WTNP1(I)/(K2D2*WNP1(1,I)) )*Fv2,0.0)

   !Part 9:
    Ft2=1.2*exp(-0.5*Chi*Chi)
    Prod(I) = A(I) * Cb1 * (1-Ft2) * WTNP1(I) * Shat

   !Part 10:
    R  = MR *  WTNP1(I) / (WNP1(1,I)*Shat*K2D2)
	if(R>10.) R=10.

	G  = R+Cw2*(R**6. - R)
	Fw = G*( (1+Cw36) / (G**6. + Cw36) )**(1/6.)

   !Part 11:
    Dest(i) = MR * A(I) * (Cw1*Fw-(Cb1/Kei*Kei)*Ft2)  *  WTNP1(I)*WTNP1(I) / (WNP1(1,I)*D2)
 End do
!*********************************************************************************************
 End
!###########################################################################################

