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
!// Chief Developer: M. Namvar, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. A. Jahangirian, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: M. Namvar, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: S. Kavoosi, Mechanical Eng., Amirkabir University of Technology          //!
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Solver_DANSInviscid(DIM,NP,NC,NF,NR,NFR,BC,IDS,X,Y,Minf,ALF,ERmx,CFLx,NRKS,NWrite,RKJ,&
							       NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
							       Xc,Yc,NX,NY,DA,A,GM,R0,P0,C0,U0,V0,&
							       WNP1,P,WB) 
 Implicit None
!**********************************************************************************************************
INTEGER                      ,INTENT(IN)    ::DIM
INTEGER                      ,INTENT(IN)    ::NP
INTEGER                      ,INTENT(IN)    ::NC
INTEGER                      ,INTENT(IN)    ::NF
INTEGER                      ,INTENT(IN)    ::NR
INTEGER,DIMENSION(1:100)     ,INTENT(IN)    ::NFR
INTEGER,DIMENSION(1:100)     ,INTENT(IN)    ::BC
INTEGER,DIMENSION(1:4,1:DIM) ,INTENT(IN)    ::IDS
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::X
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::Y
Real(8)                      ,INTENT(IN)    ::MINF
Real(8)                      ,INTENT(IN)    ::ALF
Real(8)                      ,INTENT(IN)    ::ERMX
Real(8)                      ,INTENT(IN)    ::CFLX
INTEGER                      ,INTENT(IN)    ::NRKS
INTEGER                      ,INTENT(IN)    ::NWRITE
Real(8) ,DIMENSION(1:5)      ,INTENT(IN)    ::RKJ
INTEGER                      ,INTENT(IN)    ::NF1
INTEGER                      ,INTENT(IN)    ::NF2
INTEGER                      ,INTENT(IN)    ::NFW1
INTEGER                      ,INTENT(IN)    ::NFW2
INTEGER                      ,INTENT(IN)    ::NFF1
INTEGER                      ,INTENT(IN)    ::NFF2
INTEGER                      ,INTENT(IN)    ::NFI1
INTEGER                      ,INTENT(IN)    ::NFI2
INTEGER                      ,INTENT(IN)    ::NFS1
INTEGER                      ,INTENT(IN)    ::NFS2
INTEGER                      ,INTENT(IN)    ::NFO1
INTEGER                      ,INTENT(IN)    ::NFO2
INTEGER                      ,INTENT(IN)    ::NFIF1
INTEGER                      ,INTENT(IN)    ::NFIF2
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::XC
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::YC
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::NX
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::NY
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::DA
Real(8),DIMENSION(1:DIM)     ,INTENT(IN)    ::A
Real(8)                      ,INTENT(IN)    ::GM
Real(8)                      ,INTENT(IN)    ::R0
Real(8)                      ,INTENT(IN)    ::P0
Real(8)                      ,INTENT(IN)    ::C0
Real(8)                      ,INTENT(IN)    ::U0
Real(8)                      ,INTENT(IN)    ::V0
Real(8),DIMENSION(1:4,1:DIM) ,INTENT(INOUT) ::WNP1
Real(8),DIMENSION(1:DIM)     ,INTENT(INOUT) ::P
Real(8),DIMENSION(1:5,1:DIM) ,INTENT(OUT)   ::WB

INTEGER                                     ::I,J,Init,Ncyc,NS 
Real(8)                                     ::Co,Rm,U,V,RKco,DTmin,time_end,time_start,TT

Real(8),DIMENSION(1:4,1:DIM)                ::WC,Con
Real(8),DIMENSION(1:DIM)                    ::DT
!************************************************** Main **************************************************
!Part 1:
 Call BC_Wall(DIM,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(DIM,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(DIM,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(DIM,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(DIM,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 2:
 Ncyc = 0
 Rm   = 10.0
 TT=1.0

!Part 3:
 Do While(Rm > ERmx)

   !Part 4:
    Ncyc=Ncyc+1
     
   !Part 5:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   !Part 6:
	Call TimSTP_Inviscid(DIM,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)

   !Part 7:
    Do NS=1,NRKS
   
      !Part 8:
	   RKco=RKJ(NS)

      !Part 9:
	   Call ConMeanFlow_AUSM(DIM,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

      !Part 10:
       Do J=1,NC
		  Co = RKco*DT(J)/A(J)
          WNP1(1,J) = WC(1,J) - Co* Con(1,J) 
		  WNP1(2,J) = WC(2,J) - Co* Con(2,J)
          WNP1(3,J) = WC(3,J) - Co* Con(3,J)
          WNP1(4,J) = WC(4,J) - Co* Con(4,J)

         !Part 11:
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

       End Do
       
      !Part 12: 
       Call BC_Wall(DIM,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(DIM,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(DIM,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(DIM,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(DIM,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

     End Do !Ns	
     
	!Part 13:
 	 Call ResMass(DIM,NC,WNP1,WC,DT,Ncyc,Rm)
 
     !Part 14:
     !IF( Mod(Ncyc,NWrite)==0 )Print*,'Writing Results... ',Ncyc,Rm

 End Do !Do While

!**********************************************************************************************************
 End 
!##########################################################################################################
