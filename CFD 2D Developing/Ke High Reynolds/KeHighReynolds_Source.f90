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
!// Chief Developer: N. msnkre,Aerospace eng. Amirkabir University of Technology           //!
!// Supervisor: Dr. h. hdhrnuidn,Aerospace eng. Amirkabir University of Technology         //!
!// Date: May.,04,2018                                                                     //!
!// Developed by: N. msnkre,Aerospace Eng.,Amirkabir University of Technology              //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied,Modified and Redistributed for Non-Commercial Use.                    //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeHighReynolds_Source(Dim,NC,A,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,St)
 Implicit None
!*******************************************************************************************
INTEGER(4)                   ,INTENT(IN)  ::DIM
INTEGER(4)                   ,INTENT(IN)  ::NC
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::A
REAL(8)                      ,INTENT(IN)  ::MR
REAL(8)                      ,INTENT(IN)  ::CE1
REAL(8)                      ,INTENT(IN)  ::CE2
REAL(8),DIMENSION(1:4,1:DIM) ,INTENT(IN)  ::WNP1
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(IN)  ::WTNP1
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::MU
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::MUT
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DDUX
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DDUY
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DDVX
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DDVY
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(OUT) ::ST

INTEGER                                   ::I
REAL(8)                                   ::K,Epsilon,Rho,Pk,Pe,Txx,Txy,Tyy
!******************************************************************************************* 
!Part 1:
 Do I=1,NC

   !Part 2:
    Rho     = WNP1(1,I)
    k       = WTNP1(1,I)/Rho
    Epsilon = WTNP1(2,I)/Rho

   !Part 3:
    Txx = MR * Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) - Rho*K/1.5
    Txy = MR * Mut(I)*( DDUY(I)+DDVX(I)  )
    Tyy = MR * Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) - Rho*K/1.5

   !Part 4:
    Pe =  txx*DDUX(I) + txy*(DDUY(I)+DDVX(I)) + tyy*DDVY(I) 
    Pk = min (Pe,10.0*Rho*Epsilon)

    !Part 5:
    St(1,I) = A(I) * ( -Pk + Rho*Epsilon ) 
    St(2,I) = A(I) * ( -Ce1*Pk*Epsilon/k + Ce2*Rho*Epsilon*Epsilon/k )   

 END DO
!*******************************************************************************************
 End
!###########################################################################################

