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
 Subroutine KeHighReynolds_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,Rokinf,Roeinf,WB,WTNP1,IDS,WTB)
 Implicit None
!*******************************************************************************************
INTEGER                      ,INTENT(IN)   ::DIM
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)   ::NX
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)   ::NY
INTEGER                      ,INTENT(IN)   ::NF
INTEGER                      ,INTENT(IN)   ::NFS1
INTEGER                      ,INTENT(IN)   ::NFS2
INTEGER                      ,INTENT(IN)   ::NFO1
INTEGER                      ,INTENT(IN)   ::NFO2
INTEGER                      ,INTENT(IN)   ::NFW1
INTEGER                      ,INTENT(IN)   ::NFW2
INTEGER                      ,INTENT(IN)   ::NFI1
INTEGER                      ,INTENT(IN)   ::NFI2
INTEGER                      ,INTENT(IN)   ::NFF1
INTEGER                      ,INTENT(IN)   ::NFF2
REAL(8)                      ,INTENT(IN)   ::ROKINF
REAL(8)                      ,INTENT(IN)   ::ROEINF
REAL(8),DIMENSION(1:5,1:DIM) ,INTENT(IN)   ::WB
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(IN)   ::WTNP1
INTEGER,DIMENSION(1:4,1:DIM) ,INTENT(IN)   ::IDS
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(OUT)  ::WTB

REAL(8)                                    ::U,V,Q
INTEGER                                    ::I,ME,P1,P2
!*******************************************************************************************
!Part 1:
 DO I=NFI1+1,NFI2
 	WTB(1,I) = Rokinf
    WTB(2,I) = Roeinf
 END do
  
!Part 2:  
 DO I=NFO1+1,NFO2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 3:
 DO I=NFW1+1,NFW2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 4:
 DO I=NFS1+1,NFS2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 5:
 DO I=NFF1+1,NFF2
     
   !Part 6:
    ME  = IDS(1,I)
    
   !Part 7:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    
   !Part 8:
    Q  = U*NX(I)+V*NY(I)   
    
   !Part 9:
    IF (Q<=0.0) then
 	WTB(1,I) = Rokinf
    WTB(2,I) = roeinf
    
    ELSE
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME) 
    
    END IF

 END do
!*******************************************************************************************
 End
!###########################################################################################
