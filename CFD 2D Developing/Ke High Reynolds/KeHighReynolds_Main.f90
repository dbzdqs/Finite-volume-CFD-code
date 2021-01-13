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
 Subroutine KeHighReynolds_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
           NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,DA,U0,V0,Mut0,MR,NRKS,Wb,WNP1,Mu,Rokinf,Roeinf,RKJ,TauWall,WTNP1,Mut)
 Implicit None
!*******************************************************************************************
INTEGER                      ,INTENT(IN)  ::DIM
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::X
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::Y
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::NX
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::NY
INTEGER                      ,INTENT(IN)  ::NC
INTEGER                      ,INTENT(IN)  ::NF
INTEGER                      ,INTENT(IN)  ::NF1
INTEGER                      ,INTENT(IN)  ::NF2
INTEGER                      ,INTENT(IN)  ::NFW1
INTEGER                      ,INTENT(IN)  ::NFW2
INTEGER                      ,INTENT(IN)  ::NFI1
INTEGER                      ,INTENT(IN)  ::NFI2
INTEGER                      ,INTENT(IN)  ::NFO1
INTEGER                      ,INTENT(IN)  ::NFO2
INTEGER                      ,INTENT(IN)  ::NFS1
INTEGER                      ,INTENT(IN)  ::NFS2
INTEGER                      ,INTENT(IN)  ::NFF1
INTEGER                      ,INTENT(IN)  ::NFF2
INTEGER                      ,INTENT(IN)  ::NP
INTEGER,DIMENSION(1:4,1:DIM),INTENT(IN)  ::IDS
INTEGER,DIMENSION(1:DIM)     ,INTENT(IN)  ::INW
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::XC
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::YC
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DW
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DT
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::A
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::DA
REAL(8)                      ,INTENT(IN)  ::U0
REAL(8)                      ,INTENT(IN)  ::V0
REAL(8)                      ,INTENT(IN)  ::MUT0
REAL(8)                      ,INTENT(IN)  ::MR
INTEGER                      ,INTENT(IN)  ::NRKS
REAL(8),DIMENSION(1:5,1:DIM) ,INTENT(IN)  ::WB
REAL(8),DIMENSION(1:4,1:DIM) ,INTENT(IN)  ::WNP1
REAL(8),DIMENSION(1:DIM)     ,INTENT(IN)  ::MU
REAL(8)                      ,INTENT(IN)  ::ROKINF
REAL(8)                      ,INTENT(IN)  ::ROEINF
REAL(8),DIMENSION(1:5)       ,INTENT(IN)  ::RKJ
REAL(8),DIMENSION(1:DIM)     ,INTENT(OUT) ::TAUWALL
REAL(8),DIMENSION(1:2,1:DIM) ,INTENT(OUT) ::WTNP1
REAL(8),DIMENSION(1:DIM)     ,INTENT(OUT) ::MUT

INTEGER                                   ::I,II,J,ME,NS
REAL(8)                                   ::RKco,K,Epsilon,Rho,Co,Sigk,Sige,Ce1,Ce2,Cmu,Tu
REAL(8),DIMENSION(1:2,1:Dim)              ::WTB,WTN,Cont,Dift,St
REAL(8),DIMENSION(1:Dim)                  ::DUY,DUX_C,DUY_C,DVX_C,DVY_C, DKX_F,DKY_F,DEPSX_F,DEPSY_F
!*******************************************************************************************
!Part 1:
 Sigk=1.0d0
 Sige=1.3
 Ce1=1.35d0
 Ce2=1.8d0
 Cmu=0.09d0
 
!Part 2:
 Do I=1,NC
    WTN(1,I) =WTNP1(1,I)
    WTN(2,I) =WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
    
   !Part 4:
	RKco=RKJ(NS)

   !Part 5:
    Call KeHighReynolds_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,Rokinf,Roeinf,WB,WTNP1,IDS,WTB)
 
    !Part 6:
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
    
   !Part 7:
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX_F,DKY_F,DEPSX_F,DEPSY_F)
    
   !Part 8:
    Do I=NFS1+1,NFS2
       DKY_F(I)=0.0
	   DEPSY_F(I)=0.0
    End do
   
   !Part 9:
	Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,WTNP1,Wb,WTB,Cont)

   !Part 10:
	Call KFi_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,MR,Mu,Mut,Dift)

   !Part 11:
    Call KeHighReynolds_Source(Dim,NC,A,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUX_C,DUY_C,DVX_C,DVY_C,St)

   !Part 12:
    Do J=1,NC

      !Part 13:
	   Co = RKco*DT(J)/A(J)
      
      !Part 14:
       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) + St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) + St(2,J) ) 
       
       !Part 15:
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)
       
      !Part 16:
       Rho     = WNP1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho
       Mut(J) = Cmu*Rho*k*k / Epsilon /MR

    End Do !J
    
   !Part 17:
    Call KeStandardWallFU(Dim,IDS,NFW1,NFW2,MR,NX,DA,NY,Mu,DW,WNP1,WTNP1,TauWall)

   !Part 18:
    Do J=NFW1+1,NFW2
        
      !Part 19: 
       ME      = IDS(1,J)
       Rho     = WNP1(1,ME)
       K       = WTNP1(1,ME)/Rho
       Epsilon = WTNP1(2,ME)/Rho
       
      !Part 20:
       Mut(ME) = Cmu*Rho*k*k / Epsilon /MR
       
      !Part 21:
       if(Mut(I)<0.0) Mut(i)=Mut0
     
    End Do !J

 End Do !NS
!*******************************************************************************************
 End
!###########################################################################################
