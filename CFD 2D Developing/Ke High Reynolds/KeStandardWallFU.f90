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
 Subroutine KeStandardWallFU(Dim,IDS,NFW1,NFW2,MR,NX,DA,NY,Mu,DW,WNP1,WTNP1,TauWall)
 Implicit None
!*******************************************************************************************
INTEGER                      ,INTENT(IN)    :: DIM
INTEGER,DIMENSION(1:4,1:DIM) ,INTENT(IN)    :: IDS
INTEGER                      ,INTENT(IN)    :: NFW1
INTEGER                      ,INTENT(IN)    :: NFW2
REAL(8)                      ,INTENT(IN)    :: MR
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)    :: NX
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)    :: DA
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)    :: NY
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)    :: MU
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)    :: DW
REAL(8),DIMENSION(1:4,1:DIM) ,INTENT(IN)    :: WNP1
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(INOUT) :: WTNP1
REAL(8),DIMENSION(1:DIM)     ,INTENT(OUT)   :: TAUWALL

REAL(8)                                     ::Tx,Ty,Rho,U,V,Vis,Yn,Ut,absUt,KUt,EYp,UtauN,Yp,utau
INTEGER                                     ::I,J,ME
!******************************************************************************************* 
!Part 1:
 Do I=NFW1+1,NFW2
     
   !Part 2:
    ME = IDS(1,I)
   
   !Part 3:
    Tx = -NY(I)/DA(I)
    Ty =  NX(I)/DA(I)

   !Part 4:
    Rho = WNP1(1,ME)
    U   = WNP1(2,ME)/Rho
    V   = WNP1(3,ME)/Rho

   !Part 5: 
    Vis = Mu(ME)
    Yn  = DW(ME)

   !Part 6:
    Ut  = U*Tx + V*Ty
    absUt = abs(Ut)

   !Part 7:
    utau = dsqrt( MR* (Vis*absUt) / (Rho*Yn) )

   !Part 8:
    KUt = 0.41*absUt
    EYp = 9.0*Yn*Rho/Vis/MR
    
   !Part 9:
    Do J=1,20
       IF(Utau < 0.000001) Utau = 0.000001
       UtauN = KUt/Dlog(EYp*Utau)
       Utau = UtauN
    End Do

   !Part 10:
    Yp = Rho*Utau*Yn/Vis/MR
 
   !Part 11:
    IF(Yp<12.) utau = dsqrt( MR* (Vis*absUt) / (Rho*Yn) )
 
   !Part 12:
    TauWall(I) = Rho*Utau*Utau*Ut/absUt /MR 

   !Part 13:
    WTNP1(1,ME) = Rho*utau*utau/0.3
    WTNP1(2,ME) = Rho*utau*utau*utau/(0.42*Yn)

 END DO
 
!*******************************************************************************************
 End
!###########################################################################################

