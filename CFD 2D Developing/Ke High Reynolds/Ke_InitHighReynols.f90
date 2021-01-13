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
 Subroutine Ke_InitHighReynols(Dim,NC,U0,V0,Mu0,MR,WTNP1,Mut,WNP1,Mu,Mut0,rokinf,roeinf) 
 Implicit None
!*******************************************************************************************
INTEGER                        ,INTENT(IN)  ::DIM
INTEGER                        ,INTENT(IN)  ::NC
REAL(8)                        ,INTENT(IN)  ::U0
REAL(8)                        ,INTENT(IN)  ::V0
REAL(8)                        ,INTENT(IN)  ::MU0
REAL(8)                        ,INTENT(IN)  ::MR
REAL(8),DIMENSION(1:2,1:DIM)   ,INTENT(OUT) ::WTNP1
REAL(8),DIMENSION(1:DIM)       ,INTENT(OUT) ::MUT
REAL(8),DIMENSION(1:4,1:DIM)   ,INTENT(IN)  ::WNP1
REAL(8),DIMENSION(1:DIM)       ,INTENT(IN)  ::MU
REAL(8)                        ,INTENT(OUT) ::MUT0
REAL(8)                        ,INTENT(OUT) ::ROKINF
REAL(8)                        ,INTENT(OUT) ::ROEINF

INTEGER                                     ::J,ME
REAL(8)                                     ::Tu,Cmu
!*******************************************************************************************	
!Part 1:
 Cmu=0.09d0
 Tu=0.01
 Mut0 = Tu*Mu0
 rokinf=0.000025*(U0*U0+V0*V0) 
 roeinf=cmu*rokinf*rokinf/(Mut0*MR)

!Part 2:
 Do J=1,NC
    WTNP1(1,J) = rokinf 
    WTNP1(2,J) = roeinf
    Mut(J)     = Mut0
 End Do
!*******************************************************************************************
 End
!###########################################################################################
    
    
