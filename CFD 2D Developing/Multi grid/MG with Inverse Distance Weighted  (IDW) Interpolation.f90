!****************************************************************************
!
!  PROGRAM: MULTIGRID
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

!/////////////////////////////////////////////////////////////////////////////////////////////////////////!
!// Multigrid                                                                                  //!
!// Date :         January/07/2017                                                                      //!
!// Developed by : M. Hashemi, Iran, Tehran, Mohammadhashemi@ut.ac.ir                                   //!
!//                                                                                                     //!
!// an Inviscid 2D Flow Solver                                                                          //!
!// Features: 1- 2D                                                                                     //!
!//           2- Using Unstructured Mesh                                                                //!
!//           3- Edge Based Data Structured                                                             //!
!//           4- Cell Based Conrol Volume                                                               //!
!//           5- Inviscid Flow                                                                          //!
!//           6- Convection Terms is Discritized by ECUSP Scheme                                        //!
!//           7- Transient Term is Discritized by Runge-Kutta Explicit                                  //!
!//           8- Density Based                                                                          //!
!//                                                                                                     //!
!// The Program Is Available Through The Website: www.MarketCode.ir                                     //!
!// It May be Copied, Modified, and Redistributed for Non-Commercial Use.                               //!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////!
!**********************************************************************************************************
 program Multigrid
 Implicit None
!===============================
 Integer,Parameter::Dim=120000
!===============================
 Integer::I,J,K,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,NRKS,NWrite,&
          NR,Init,Ncyc,NS,counter,C,NSN,NMF,NE,ME,NOL,FTC,NON,MGcycle
 Real(8)::GM,Co,ALF,R0,P0,C0,U0,V0,E0,ERmx,CFLx,Minf,Rm,U,V,RKco,DTmin,ERmx1,epsilon,ER1,ER2,ER3,ER4
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS,Corn
 Integer,Dimension(1:3,1:Dim)::Corn1,Corn2,CornC,NeibC
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,RES,RESM0,RESM1,RESM2,Error,SolM1,WC1,FFM2
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,A,NX,NY,DA,DT,P,X2,Y2,A1,X1,Y1
 Real(8),Dimension(1:5,1:Dim)::WB
 real :: start, finish
!************************************************** Main **************************************************
call cpu_time(start)
!Part 1:
 NMF=1
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)

!Part 2:
 Call Read_SettingV1(Minf,ALF,ERmx,CFLx,NRKS,NWrite,Init)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 4:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 5:
 Call InitMeanFlow_Inviscid(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,WNP1)
 
!Part 6: 
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

!Part 7:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 
!Part 8:
 Do I=1,2
       
   !Part 9:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   !Part 10:
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)
   
   !Part 11:
    Do NS=1,NRKS
   
      !Part 12:
	   RKco=1.0/(NRKS-NS+1)

      !Part 13:
	   Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

      !Part 14:
       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co* Con(1,J)
		  WNP1(2,J) = WC(2,J) - Co* Con(2,J)
          WNP1(3,J) = WC(3,J) - Co* Con(3,J)
          WNP1(4,J) = WC(4,J) - Co* Con(4,J)

	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

       End Do

      !Part 15: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)

    End Do !Ns	
    
 End Do !Do

 !Part 16:
 Do I=1,NC  
     SolM1(1,I)=WNP1(1,I)
     SolM1(2,I)=WNP1(2,I)
     SolM1(3,I)=WNP1(3,I)
     SolM1(4,I)=WNP1(4,I)
 End DO

!Part 17: 
Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

!Part 18:
Do I=1,NC
    RES(1,I)=Con(1,I)
    RES(2,I)=Con(2,I)
    RES(3,I)=Con(3,I)
    RES(4,I)=Con(4,I)
End Do

!Part 19:
NOL=1
FTC=1
Call Interpolation2D(Dim,WNP1,RES,Error,NOL,FTC)

 !Part 20:
 NMF=2
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
 
 !Part 21:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
 !Part 22:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 
 !Part 23:
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

 !Part 24:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 
 !Part 25:
 Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

 !Part 26:
 Do I=1,NC
    RESM2(1,I)=Con(1,I)
    RESM2(2,I)=Con(2,I)
    RESM2(3,I)=Con(3,I)
    RESM2(4,I)=Con(4,I)  
 End Do
 
 !Part 27:
 Do I=1,NC
    FFM2(1,I)=RES(1,I)-RESM2(1,I)
    FFM2(2,I)=RES(2,I)-RESM2(2,I)
    FFM2(3,I)=RES(3,I)-RESM2(3,I)
    FFM2(4,I)=RES(4,I)-RESM2(4,I)
 End Do

 !Part 28:
 Do I=1,NC
    WC1(1,I)=WNP1(1,I)
    WC1(2,I)=WNP1(2,I)
    WC1(3,I)=WNP1(3,I)
    WC1(4,I)=WNP1(4,I)
 End Do

!Part 29:
Rm=10
Ncyc=0
ERmx1=-6.5

!Part 30:
Do While(Rm > ERmx1)

   !Part 31:
    Ncyc=Ncyc+1
     
   !Part 32:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J) 
    End Do

   !Part 33:
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)
    
   !Part 34:
    Do NS=1,NRKS
   
      !Part 35:
	   RKco=1.0/(NRKS-NS+1)

      !Part 36:
	   Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

      !Part 37:
       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co* (Con(1,J))
		  WNP1(2,J) = WC(2,J) - Co* (Con(2,J))
          WNP1(3,J) = WC(3,J) - Co* (Con(3,J))
          WNP1(4,J) = WC(4,J) - Co* (Con(4,J))

          
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

       End Do

      !Part 38: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)

    End Do !Ns	

    !Part 39:
 	 Call ResMass1(Dim,NC,WNP1,WC,DT,Rm)
     Print*,Ncyc,Rm

End Do !Do While


!Part 40:
Do I=1,NC
    Error(1,I)=WNP1(1,I)-WC1(1,I)
    Error(2,I)=WNP1(2,I)-WC1(2,I)
    Error(3,I)=WNP1(3,I)-WC1(3,I)
    Error(4,I)=WNP1(4,I)-WC1(4,I)
End Do

 !Part 41:
 NOL=1
 FTC=0
 Call Interpolation2D(Dim,WNP1,RES,Error,NOL,FTC)
 
 !Part 42:
 NMF=1
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
 
 !Part 43:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
 !Part 44:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 
 
!Part 45:
Do I=1,NC
    WNP1(1,I)=SolM1(1,I)+Error(1,I)
    WNP1(2,I)=SolM1(2,I)+Error(2,I)
    WNP1(3,I)=SolM1(3,I)+Error(3,I)
    WNP1(4,I)=SolM1(4,I)+Error(4,I)
End Do


 !Part 46:
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

 !Part 47:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)

!Part 48:
Rm=10
Ncyc=0
Do While(Rm > ERmx)
    Ncyc=Ncyc+1
   
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)

    Do NS=1,NRKS
   
	   RKco=1.0/(NRKS-NS+1)

	   Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)

       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co* Con(1,J)
		  WNP1(2,J) = WC(2,J) - Co* Con(2,J)
          WNP1(3,J) = WC(3,J) - Co* Con(3,J)
          WNP1(4,J) = WC(4,J) - Co* Con(4,J)

         
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

       End Do

       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)

    End Do !Ns	
    
    Call ResMass(Dim,NC,WNP1,WC,DT,Rm)
     Print*,Ncyc,Rm

     !Part 49:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
	  Call Write_ResultsV1(Dim,NFW1,NFW2,NC,NP,NF,X,Y,IDS,GM,Minf,WNP1,P)
	 End If
End Do !while


 Call Write_ResultsV1(Dim,NFW1,NFW2,NC,NP,NF,X,Y,IDS,GM,Minf,WNP1,P)
 call cpu_time(finish)
 print *, finish-start  
 
 pause
!**********************************************************************************************************
 End program Multigrid
!##########################################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Interpolation between two Meshes by intersection of triangles           //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: January, 07, 2017                                                              //!
!// Developed by: M. Hashemi, Iran, Tehran, Mohammadhashemi@ut.ac.ir                     //!
!// Doc ID: MC2F100F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
Subroutine Interpolation2D(Dim,WNP1,RES,Error,NOL,FTC)
implicit none
!***********************************************************************************************************
Intent (IN)::Dim,NOL
Intent (InOut)::WNP1,RES,Error
!===============================
 Integer::Dim,I,J,K,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NS,counter,counter1,NMF,NP1,NC1,NF_1,NF_2,NR1,NP2,NC2,NR2,P1,P2,P3,ME,NE,h1,h2,counter10,counter5,counter6,num1,counter3,counter9,counter11,counter13,NP2p,NP1p,NOL,FTC,PP
 Real(8)::F1,F2,F3,F_p1,F_p2,F_p3,F_s1,F_s2,F_s3,F_h,F_h1,F_h2,F_h3,collisionX,collisionY,M1,M2,Multi1,Multi2,a12,a13,a23,d,d_total,subX,subY
 Integer,Dimension(1:100)::NFR,BC,NFR1,BC1,NFR2,BC2
 Integer,Dimension(1:Dim)::SPN
 Integer,Dimension(1:3)::P_1,P_2
 Integer,Dimension(1:4,1:Dim)::IDS1,IDS2,IDS,corn
 Integer,Dimension(1:3,1:Dim)::Corn1,Corn2
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,RES,Error,delta_s0,delta_s1,delta_s2,delta_sm,delta_sn,delta_sp,FUNC1,FUNC2,FUNC3,PFUNC1,PFUNC2,PFUNC3,FUNC_1,FUNC_2,FUNC_3
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,XC1,YC1,XC2,YC2,A,NX,NY,DA,DT,P,X2,Y2,A1,X1,Y1,A_t,A_h
 Real(8),Dimension(1:5,1:Dim)::WB
!************************************************************************************************************
 Do I=1,Dim
    
    FUNC1(1,I)=0;FUNC1(2,I)=0;FUNC1(3,I)=0;FUNC1(4,I)=0
    FUNC_1(1,I)=0;FUNC_1(2,I)=0;FUNC_1(3,I)=0;FUNC_1(4,I)=0
    
    FUNC2(1,I)=0;FUNC2(2,I)=0;FUNC2(3,I)=0;FUNC2(4,I)=0    
    FUNC_2(1,I)=0;FUNC_2(2,I)=0;FUNC_2(3,I)=0;FUNC_2(4,I)=0
        
    FUNC3(1,I)=0;FUNC3(2,I)=0;FUNC3(3,I)=0;FUNC3(4,I)=0    
    FUNC_3(1,I)=0;FUNC_3(2,I)=0;FUNC_3(3,I)=0;FUNC_3(4,I)=0
    
End Do
 

!Part 3:
if(NOL==1 .AND. FTC==1)then
NMF=1
else if(NOL==1 .AND. FTC==0)then
NMF=2
end if

!Part 6:
Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
NC1=NC
!Part 7:
Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 9:
Do J=1,NC1  
    Xc1(J)=Xc(J);Yc1(J)=Yc(J)
End Do


Do I=1,NC1
    FUNC1(1,I)=WNP1(1,I);FUNC1(2,I)=WNP1(2,I);FUNC1(3,I)=WNP1(3,I);FUNC1(4,I)=WNP1(4,I)
    FUNC2(1,I)=RES(1,I);FUNC2(2,I)=RES(2,I);FUNC2(3,I)=RES(3,I);FUNC2(4,I)=RES(4,I)
    FUNC3(1,I)=Error(1,I);FUNC3(2,I)=Error(2,I);FUNC3(3,I)=Error(3,I);FUNC3(4,I)=Error(4,I)
end Do

Do I=1,NC1
    PFUNC1(1,I)=FUNC1(1,I);PFUNC1(2,I)=FUNC1(2,I);PFUNC1(3,I)=FUNC1(3,I);PFUNC1(4,I)=FUNC1(4,I)
    PFUNC2(1,I)=FUNC2(1,I);PFUNC2(2,I)=FUNC2(2,I);PFUNC2(3,I)=FUNC2(3,I);PFUNC2(4,I)=FUNC2(4,I)  
    PFUNC3(1,I)=FUNC3(1,I);PFUNC3(2,I)=FUNC3(2,I);PFUNC3(3,I)=FUNC3(3,I);PFUNC3(4,I)=FUNC3(4,I)
End Do


if(NOL==1 .AND. FTC==1)then
NMF=2
else if(NOL==1 .AND. FTC==0)then
NMF=1
end if

!Part 6:
Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
NC2=NC
!Part 7:
Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 
Do J=1,NC2
    Xc2(J)=Xc(J);Yc2(J)=Yc(J)
End Do
      

Do I=1,NC2
    
    d_total=0
    
    Do J=1,NC1
        
        subX=Xc2(I)-Xc1(J)
        subY=Yc2(I)-Yc1(J)
        d=sqrt(subX*subX+subY*subY)
        
        if(d==0)then
        cycle
        end if
        
        d=d**8
        d_total=d_total+1/d
        
        FUNC_1(1,I)=FUNC_1(1,I)+FUNC1(1,J)/d
        FUNC_1(2,I)=FUNC_1(2,I)+FUNC1(2,J)/d
        FUNC_1(3,I)=FUNC_1(3,I)+FUNC1(3,J)/d
        FUNC_1(4,I)=FUNC_1(4,I)+FUNC1(4,J)/d
        
        FUNC_2(1,I)=FUNC_2(1,I)+FUNC2(1,J)/d
        FUNC_2(2,I)=FUNC_2(2,I)+FUNC2(2,J)/d
        FUNC_2(3,I)=FUNC_2(3,I)+FUNC2(3,J)/d
        FUNC_2(4,I)=FUNC_2(4,I)+FUNC2(4,J)/d
        
        FUNC_3(1,I)=FUNC_3(1,I)+FUNC3(1,J)/d
        FUNC_3(2,I)=FUNC_3(2,I)+FUNC3(2,J)/d
        FUNC_3(3,I)=FUNC_3(3,I)+FUNC3(3,J)/d
        FUNC_3(4,I)=FUNC_3(4,I)+FUNC3(4,J)/d
        
    End Do
    
        FUNC_1(1,I)=FUNC_1(1,I)/d_total
        FUNC_1(2,I)=FUNC_1(2,I)/d_total
        FUNC_1(3,I)=FUNC_1(3,I)/d_total
        FUNC_1(4,I)=FUNC_1(4,I)/d_total
        
        FUNC_2(1,I)=FUNC_2(1,I)/d_total
        FUNC_2(2,I)=FUNC_2(2,I)/d_total
        FUNC_2(3,I)=FUNC_2(3,I)/d_total
        FUNC_2(4,I)=FUNC_2(4,I)/d_total
        
        FUNC_3(1,I)=FUNC_3(1,I)/d_total
        FUNC_3(2,I)=FUNC_3(2,I)/d_total
        FUNC_3(3,I)=FUNC_3(3,I)/d_total
        FUNC_3(4,I)=FUNC_3(4,I)/d_total 
End Do

Do I=1,NC2
    FUNC1(1,I)=FUNC_1(1,I)
    FUNC1(2,I)=FUNC_1(2,I)
    FUNC1(3,I)=FUNC_1(3,I)
    FUNC1(4,I)=FUNC_1(4,I)
    
    FUNC2(1,I)=FUNC_2(1,I)
    FUNC2(2,I)=FUNC_2(2,I)
    FUNC2(3,I)=FUNC_2(3,I)
    FUNC2(4,I)=FUNC_2(4,I)
    
    FUNC3(1,I)=FUNC_3(1,I)
    FUNC3(2,I)=FUNC_3(2,I)
    FUNC3(3,I)=FUNC_3(3,I)
    FUNC3(4,I)=FUNC_3(4,I)
END Do


!Part 58:
Do I=1,NC2
    counter9=0
    Do J=1,NC1
        
    if(Xc2(I)==Xc1(J) .AND. Yc2(I)==Yc1(J))then
    counter9=counter9+1    
    end if
    
    
    if(counter9/=0)then
        
            FUNC1(1,I)=PFUNC1(1,J)
            FUNC1(2,I)=PFUNC1(2,J)
            FUNC1(3,I)=PFUNC1(3,J)
            FUNC1(4,I)=PFUNC1(4,J)
            
            FUNC2(1,I)=PFUNC2(1,J)
            FUNC2(2,I)=PFUNC2(2,J)
            FUNC2(3,I)=PFUNC2(3,J)
            FUNC2(4,I)=PFUNC2(4,J)
            
            FUNC3(1,I)=PFUNC3(1,J)
            FUNC3(2,I)=PFUNC3(2,J)
            FUNC3(3,I)=PFUNC3(3,J)
            FUNC3(4,I)=PFUNC3(4,J)
           
            Exit    
    end if
        
    End Do
End Do


Do I=1,NC2
    WNP1(1,I)=FUNC1(1,I);WNP1(2,I)=FUNC1(2,I);WNP1(3,I)=FUNC1(3,I);WNP1(4,I)=FUNC1(4,I)
    RES(1,I)=FUNC2(1,I);RES(2,I)=FUNC2(2,I);RES(3,I)=FUNC2(3,I);RES(4,I)=FUNC2(4,I)
    Error(1,I)=FUNC3(1,I);Error(2,I)=FUNC3(2,I);Error(3,I)=FUNC3(3,I);Error(4,I)=FUNC3(4,I)
end Do

end
!##########################################################################################################
!##########################################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Apply InFlow Boundary Condition                                         //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F017F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine	BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NFI1,NFI2,GM,U0,V0,P0,R0,IDS,Wnp1,NX,NY,DA,P,ALF,Minf
 Intent (Out  )::Wb

 Integer::Dim,I,NFI1,NFI2,ME
 Real(8)::GM,GM1,U,V,CC,MLoc,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,NXX,NYY,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0,ITT,ITP,ALF,HTE,RIB,&
		  HTB,EP1,EP2,EP3,C1,C2,MB,TB,Minf,UBB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::WB
!*******************************************************************************************
 ITT=(1.+0.2*Minf*Minf)
 ITP= (ITT**3.5)   * P0

 GM1= GM-1
 DO I=NFI1+1,NFI2

    ME  = IDS(1,I)

    NXX = NX(I)/DA(I)
    NYY = NY(I)/DA(I)

    U = WNP1(2,ME)/WNP1(1,ME)
    V = WNP1(3,ME)/WNP1(1,ME)

    CC = GM*P(ME)/WNP1(1,ME)

    MLoc = SQRT((U*U+V*V)/CC)

    R0B = R0
    U0B = U0
    V0B = V0
    P0B = P0
    C0B = SQRT(GM*P0B/R0B)

    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
    CE = SQRT(ABS(GM*PE/RE))

    HTE = (PE/RE)*(GM/GM1)+0.5*(UE*UE+VE*VE)

    QNE = UE*NXX+VE*NYY
    QTE =-UE*NYY+VE*NXX
    RIE =-QNE-2.*CE/GM1

	IF(MLoc<1.)then
	RIB=RIE
	HTB=HTE

	EP1=1+2/GM1
	EP2=2*RIB
	EP3=(GM1/2)*(RIB*RIB-2*HTB)

	C1=(-EP2+SQRT(EP2*EP2-4*EP1*EP3))/(2*EP1)
	C2=(-EP2-SQRT(EP2*EP2-4*EP1*EP3))/(2*EP1)

	C=MAX(C1,C2)

	UBB= (2*C/GM1)+RIB

    MB=UBB/C

	PB=ITP*(1+(GM1/2)*MB**2)**(-GM/GM1)

	TB=ITT*(PB/ITP)**(GM1/GM)

	RB=GM*PB/TB

    UB=Cos(ALF)*UBB
	VB=Sin(ALF)*UBB

	END IF

	REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)

	Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

 END do
!*******************************************************************************************
 End
!########################################################################################### 
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Apply Riemann Boundary Condition                                        //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F015F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NFF1,NFF2,GM,U0,V0,P0,R0,C0,IDS,Wnp1,NX,NY,DA,P
 Intent (Out  )::Wb

 Integer::Dim,I,NFF1,NFF2,ME,P1,P2
 Real(8)::GM,GM1,U,V,CC,MLoc,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,STH,CTH,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::Wb
!*******************************************************************************************	
  
 GM1= GM-1.

 DO I=NFF1+1,NFF2

    ME  = IDS(1,I)
    P1  = IDS(3,I)
    P2  = IDS(4,I)

    STH = NX(I)/DA(I)
    CTH = NY(I)/DA(I)

    U = WNP1(2,ME)/WNP1(1,ME)
    V = WNP1(3,ME)/WNP1(1,ME)

    CC = GM*P(ME)/WNP1(1,ME)

    MLoc = SQRT((U*U+V*V)/CC)

    U0B = U0
    V0B = V0
    P0B = P0
    R0B = R0
    C0B = SQRT(GM*P0B/R0B)

    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
    CE = SQRT(ABS(GM*PE/RE))

    QN0 = U0B*STH+V0B*CTH
    QT0 =-U0B*CTH+V0B*STH
    RI0 = QN0 - 2.*C0B/GM1
    S0  = P0B/R0B**GM

    QNE = UE*STH+VE*CTH
    QTE =-UE*CTH+VE*STH
    RIE = QNE+2.*CE/GM1
    SE  = PE/RE**GM

    IF(MLoc<=1.)then

     QNN = (RIE+RI0)/2.0
     C   = 0.25*GM1*(RIE-RI0)
     QTT = QT0
     SB  = S0
     IF(QNN>0.0)then
	  QTT = QTE
      SB  = SE
	 Endif

    Elseif(MLoc>1.)then

     QNN = QN0
     C   = C0
     QTT = QT0
     SB  = S0
     
	 IF(QNN>0.0)then
      QNN = QNE
      C   = CE
      QTT = QTE
      SB  = SE
	 Endif

	Endif

    RB = (C*C/GM/SB)**(1.0/GM1)
    UB = QNN*STH-QTT*CTH
    VB = QTT*STH+QNN*CTH
    PB = C*C*RB/GM

    REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)

    Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

 End do

!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Apply Subsonic OutFlow Boundary Condition on faces ajacent to Viscous Wall!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F016F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine	 BC_SubOutFlow(Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P,WB)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NFO1,NFO2,NX,NY,DA,IDS,GM,P0,WNP1,P
 Intent(Out  )::WB

 Integer::Dim,I,ME,NFO1,NFO2
 Real(8)::GM,GM1,NXX,NYY,U,V,CC,MLoc,PE,RE,UE,VE,TE,RB,UB,VB,PB,TB,REB,P0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::WB
!*******************************************************************************************	
 GM1= GM-1

 DO I=NFO1+1,NFO2

    ME  = IDS(1,I)

    U = WNP1(2,ME)/WNP1(1,ME)
    V = WNP1(3,ME)/WNP1(1,ME)

    CC = GM*P(ME)/WNP1(1,ME)
	MLoc = SQRT((U*U+V*V)/CC)

    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
	TE = PE*GM/RE

	IF(MLoc<1.)then
	PB=P0
	UB=UE
	VB=VE
	TB=TE
	RB=PB*GM/TB

	Elseif(MLoc>1.)then
	PB=PE
	UB=UE
	VB=VE
	TB=TE
	RB=PB*GM/TB

	END IF

	REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)
    
	Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

 END do
!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Apply Symmetry Boundary Condition                                       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F018F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine	 BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NFS1,NFS2,GM,U0,V0,P0,R0,C0,IDS,Wnp1,NX,NY,DA,P
 Intent (Out  )::Wb

 Integer::Dim,I,NFS1,NFS2,ME,P1,P2
 Real(8)::GM,GM1,U,V,CC,MLoc,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,STH,CTH,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,QNB,QTB,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0,ITT,ITP,AII,&
		  HTE,RIB,HTB,EP1,EP2,EP3,C1,C2,MB,TB,TE
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::Wb
!*******************************************************************************************
 GM1= GM-1
 DO I=NFS1+1,NFS2

    ME  = IDS(1,I)

    STH = NX(I)/DA(I)
    CTH = NY(I)/DA(I)

    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
	TE=PE*GM/RE
    QNE = UE*STH+VE*CTH
    QTE =-UE*CTH+VE*STH

	RB=RE
	QNB=0
	QTB=QTE
	PB=PE
	TB=TE

	UB=STH*QNB-CTH*QTB
	VB=CTH*QNB+STH*QTB
	REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)

	Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

	END do

!*******************************************************************************************
 End
!###########################################################################################
 !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Apply Wall Boundary Condition on Inviscid and Viscouse Wall Faces       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F011F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,GM,P
 Intent(Out  )::WB

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::GM,GM1,PB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:5,1:Dim)::WB
!*******************************************************************************************	
  
 GM1= GM-1.

 DO I=NFW1+1,NFW2

    ME = IDS(1,I)

    PB      = P(ME)

    WB(1,I) = 1.0
    WB(2,I) = 0.0
    WB(3,I) = 0.0
    WB(4,I) = PB/GM1
    WB(5,I) = PB

 End do
!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Convection Terms of 2D Mean Flow Equations Using AUSM     //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
!*******************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

   !Part 11:
    a_L  = Dsqrt(GM*P(L)/WNP1(1,L))
	a_R  = Dsqrt(GM*P(R)/WNP1(1,R))

	M_L = (WNP1(2,L)*NXX+WNP1(3,L)*NYY) / (WNP1(1,L) * a_L)
	M_R = (WNP1(2,R)*NXX+WNP1(3,R)*NYY) / (WNP1(1,R) * a_R)

   !Part 12:
	IF(Dabs(M_L)>1.)Then
	 M_Plus = 0.5*(M_L+Dabs(M_L))
	 P_Plus = 0.5*(M_L+Dabs(M_L))/(M_L)
	Else
     M_Plus = 0.25*(M_L+1.)*(M_L+1.)
	 P_Plus = 0.25*(M_L+1.)*(M_L+1.)*(2.-M_L)
	End If

   !Part 13:
	IF(Dabs(M_R)>1.)Then
	 M_Minus = 0.5*(M_R-Dabs(M_R))
	 P_Minus = 0.5*(M_R-Dabs(M_R))/(M_R)
	Else
     M_Minus =-0.25*(M_R-1.)*(M_R-1.)
	 P_Minus = 0.25*(M_R-1.)*(M_R-1.)*(2.+M_R)
	End If

   !Part 14:
    Mm = M_Plus+M_Minus
	Pm = P(L)*P_Plus + P(R)*P_Minus 

   !Part 15:
	If(Mm<=0.)Then
     Ro = WNP1(1,R)       * a_R
	 RU = WNP1(2,R)       * a_R
	 RV = WNP1(3,R)       * a_R
	 RH =(WNP1(4,R)+P(R)) * a_R
	Else
     Ro = WNP1(1,L)       * a_L
	 RU = WNP1(2,L)       * a_L
	 RV = WNP1(3,L)       * a_L
	 RH =(WNP1(4,L)+P(L)) * a_L
	Endif

   !Part 16:
    F1 = ( Mm * Ro          ) * DAA
    F2 = ( Mm * RU + Pm*NXX ) * DAA
    F3 = ( Mm * RV + Pm*NYY ) * DAA
    F4 = ( Mm * RH          ) * DAA

   !Part 17:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4

   !Part 18:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    
 End Do
 
!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate Geometrical Values of 2D Mesh.                                //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F026F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NF1,NF2,NF,NC,IDS,X,Y
 Intent(Out  )::Xc,Yc,NX,NY,DA,A

 Integer::Dim,I,P1,P2,P3,NF1,NF2,NC,NF,ME,NE
 Real(8)::SumX,SumY,DArea,DX,DY,DL
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,A,Xc,Yc,NX,NY,DA
!*******************************************************************************************	
!Part 1:
 DO I=1,NC
    A(I)  = 0.0
    Xc(I) = 0.0
    Yc(I) = 0.0
 End do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
	P1 = IDS(3,I)
    P2 = IDS(4,I)

   !Part 4:
    DArea = X(P1)*Y(P2) - X(P2)*Y(P1)

    SumX = X(P1) + X(P2)
    SumY = Y(P1) + Y(P2)

   !Part 5:
    A(ME) = A(ME) + DArea

    Xc(ME) = Xc(ME) + SumX*DArea
    Yc(ME) = Yc(ME) + SumY*DArea

   !Part 6:
    A(NE) = A(NE) - DArea

    Xc(NE) = Xc(NE) - SumX*DArea
    Yc(NE) = Yc(NE) - SumY*DArea
 End do

!Part 7:
 DO I=NF2+1,NF

    ME = IDS(1,I)
	P1 = IDS(3,I)
    P2 = IDS(4,I)

    DArea = X(P1)*Y(P2) - X(P2)*Y(P1)

    SumX = X(P1) + X(P2)
    SumY = Y(P1) + Y(P2)

    A(ME) = A(ME) + DArea

    Xc(ME) = Xc(ME) + SumX*DArea
    Yc(ME) = Yc(ME) + SumY*DArea

 End do

!Part 8:
 DO I=1,NC
    A(I)  = A(I)  / 2
    Xc(I) = Xc(I) / (6*A(I))
    Yc(I) = Yc(I) / (6*A(I))
 End do

!Part 9:
 DO I=1,NF

   !Part 10:
	P1 = IDS(3,I)
    P2 = IDS(4,I)

   !Part 11:
    NX(I) = Y(P2)-Y(P1)
    NY(I) = X(P1)-X(P2)
    DA(I) = Dsqrt(NX(I)*NX(I) + NY(I)*NY(I))

 End do

!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Initialize all of the Parameters Contribute in Inviscid Solver        //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F010F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: WWW.IRCSE.IR                          //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine InitMeanFlow_Inviscid(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,WNP1)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,Init,NC,Minf,ALF
 Intent(Out  )::GM,R0,P0,C0,U0,V0,WNP1

 Integer::Dim,Init,I,NC
 Real(8)::PI,ALF,Minf,GM,R0,P0,C0,U0,V0,E0,ALFA
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*******************************************************************************************	
!Part 1:
 PI  = 4.0*Atan(1.0)
 ALFA= ALF*PI/180.

!Part 2:
 GM  = 1.4

!Part 3:
 R0 = 1.0

!Part 4:
 P0 = 1.0/GM  

!Part 5:
 C0 = SQRT(GM*P0/R0)

!Part 6:
 U0 = Minf*C0*COS(ALFA)
 V0 = Minf*C0*SIN(ALFA)

!Part 7:
 E0 = P0/(R0*(GM-1))+ 0.5*(U0*U0 + V0*V0) 
 
!Part 8:
 DO I=1,NC
    WNP1(1,I) = R0
    WNP1(2,I) = R0*U0
    WNP1(3,I) = R0*V0
    WNP1(4,I) = R0*E0
 end do

!Part 9:
 IF(Init==1)Then
  Open(1,File='SolutionData.txt')
  Do I=1,NC
	 Read(1,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I)
  End Do       
  close(1)
 Endif
!*******************************************************************************************
 End
!###########################################################################################
 !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Renumbering Faces According to Boundary Conditions and                  //!
!//              Determining the Index of first and last Faces of Each Boundary Conditions/!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F022F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NR,NF
 Intent(Out  )::NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2
 Intent(InOut)::NFR,BC,IDS

 Integer::Dim,J,JJ,J1,I,SF,N,M,NR,NF,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,&
          NF1,NF2,NFIF1,NFIF2,NFN,NFW,NFF,NFI,NFO,NFS,NFIF
 Integer,Dimension(1:100)::NFR,TNFR,BC,TBC
 Integer,Dimension(1:4,1:Dim)::IDS,TIDS
!*******************************************************************************************
!Part 1:
 Do J=1,NF
	Do J1=1,4
       TIDS(J1,J) = IDS(J1,J)
    End do
 End Do

!Part 2:
 Do J=1,NR
    TNFR(J) = NFR(J)
	TBC(J)  = BC(J) 
 End do

!Part 3:
 N=0
 M=0
 Do JJ=1,10

   !Part 4:
    SF=0
    Do J=1,NR
       IF(TBC(J)==JJ)Then

	    Do I=SF+1,SF+TNFR(J)

           N=N+1
           Do J1=1,4
              IDS(J1,N) = TIDS(J1,I)
		   End do

        Enddo
		
	   !Part 5:		   
	    M=M+1
	    NFR(M) = TNFR(J)
		BC(M)  = TBC(J)  

	   Endif
	   SF=SF+TNFR(J)
    End Do

 End Do

!Part 6:
 NFN  = 0
 NFW  = 0
 NFF  = 0
 NFI  = 0
 NFO  = 0
 NFS  = 0
 NFIF = 0
 Do J=1,NR
    IF( BC(J)==1 ) NFN  = NFN  + NFR(J)
    IF( BC(J)==2 ) NFW  = NFW  + NFR(J)
    IF( BC(J)==3 ) NFF  = NFF  + NFR(J)
    IF( BC(J)==4 ) NFI  = NFI  + NFR(J)
    IF( BC(J)==5 ) NFO  = NFO  + NFR(J)
    IF( BC(J)==6 ) NFS  = NFS  + NFR(J)
    IF( BC(J)==7 ) NFIF = NFIF + NFR(J)
 End Do

!Part 7:
 NF1=0
 NF2=NF1+NFN

 NFW1=NF2
 NFW2=NFW1+NFW

 NFF1=NFW2
 NFF2=NFF1+NFF

 NFI1=NFF2
 NFI2=NFI1+NFI

 NFO1=NFI2
 NFO2=NFO1+NFO

 NFS1=NFO2
 NFS2=NFS1+NFS

 NFIF1=NFI2
 NFIF2=NFIF1+NFIF
!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Read Edge Based Mesh From 'Mesh.Txt' File                            //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: November, 6, 2015                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC5F003F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
Subroutine Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NMF
 Intent(Out  )::NP,NC,NF,NR,NFR,BC,IDS,X,Y

 Integer::Dim,I,J,J1,JJ,NP,NC,NF,NR,SFace,FaceType,MeshDim,NMF
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
!*******************************************************************************************
!Part 1:
 IF(NMF==1)Then
 Open(3,File='Mesh1.Txt')
Else IF(NMF==2)Then
 Open(3,File='Mesh2.Txt')
Else IF(NMF==3)Then
 Open(3,File='Mesh3.Txt')
End IF

!Part 2:
 Read(3,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(3,*) NP    

!Part 4:
 Read(3,*) NC

!Part 5:
 Read(3,*) NF

!Part 6:
 Read(3,*) NR

!Part 7:
 Read(3,*)
 Do J=1,NR
    Read(3,*) NFR(J) , BC(J)
 End Do 

!Part 8:
 Read(3,*)
 Do J=1,NF  
    Read(3,*) FaceType,IDS(1,J),IDS(2,J),IDS(3,J),IDS(4,J)
 End Do

!Part 9:
 Read(3,*)        
 Do J=1,NP
    Read(3,*) X(J),Y(J)
 End Do

 Close(3)
!*******************************************************************************************
End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Take all of the Necessary setting Parameters from User                //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F001F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: WWW.IRCSE.IR                          //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Read_SettingV1(Minf,ALF,ERmx,CFLx,NRKS,NWrite,Init)
 Implicit None
!******************************************************************************************* 
 Intent(Out  )::Minf,ALF,ERmx,CFLx,NRKS,NWrite,Init

 Real(8)::Minf,ALF,ERmx,CFLx
 Integer::NRKS,NWrite,Init
!*******************************************************************************************	
 Open(1,File='Setting.Txt')

!Part 1:
 Read(1,*) Minf

!Part 2:
 Read(1,*) ALF

!Part 3:
 Read(1,*) ERmx

!Part 4:
 Read(1,*) CFLx

!Part 5:
 Read(1,*) NRKS

!Part 6: 
 Read(1,*) NWrite

!Part 7: 
 Read(1,*) Init
 
 Close(1)
!*******************************************************************************************
 End
!###########################################################################################
 !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Calculate Residual of Mass Equation                                   //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F004F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine ResMass(Dim,NC,WNP1,WN,DT,Rm)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,WNP1,WN ,DT
 Intent(Out  )::Rm

 Integer::Dim,I,NC
 Real(8)::Rm
 Real(8),Dimension(1:4,1:Dim)::Wnp1,WN
 Real(8),Dimension(1:Dim)::DT
!*******************************************************************************************	
 Open(11,File='ResMass.Plt')

!Part 1:
 Rm = 0.0

!Part 2:
 DO I=1,NC
    Rm = Rm + Abs( Wnp1(1,I)-WN(1,I) ) / WN(1,I)
 End Do

!Part 3:
 Rm = Rm / NC

!Part 4:
 Rm = Dlog10(Rm+1.0E-16)

 Write(11,*) Rm
!*******************************************************************************************
 End
!###########################################################################################
 !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Calculate Time Step of 2D Inviscid Flow (Explicit Scheme)             //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: NC5F012F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: WWW.IRCSE.IR                          //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,CFLx,IDS,A,P,GM,WNP1,Wb
 Intent(Out  )::DT

 Integer::Dim,I,NC,ME,NE,NF1,NF2,NF
 Real(8)::U,V,T1,R,R1,R2,C,DX,CFLx,GM
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::DT,A,P,NX,NY,DA
!*******************************************************************************************	
!Part 1:
 Do I=1,NC
    DT(I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    R = WB(1,I)
    U = WB(2,I)/R
    V = WB(3,I)/R

   !Part 5:
    C = SQRT( ABS( GM*WB(5,I)/R ) )

   !Part 6:
    DT(ME) = DT(ME) + ABS(U*NX(I) + V*NY(I)) + C*DA(I)

 End Do

!Part 7:
 DO I=NF1+1,NF2
 
   !Part 8:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 9:
    R1 = WNP1(1,ME)
	R2 = WNP1(1,NE)
    U = 0.5*( WNP1(2,ME)/R1 + WNP1(2,NE)/R2 )
    V = 0.5*( WNP1(3,ME)/R1 + WNP1(3,NE)/R2 )

   !Part 10:
    C = SQRT( ABS( GM*(P(ME)+P(NE))/(R1+R2) ) )

   !Part 11:
    T1 = ABS(U*NX(I) + V*NY(I)) + C*DA(I)

   !Part 12:
    DT(ME) = DT(ME) + T1
    DT(NE) = DT(NE) + T1
 End Do

!Part 13:
 DO I=1,NC
    DT(I) = CFLx*A(I)/DT(I)
 End Do
!*******************************************************************************************
 End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Write the Results                                                    //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F002F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Write_ResultsV1(Dim,NFW1,NFW2,NC,NP,NF,X,Y,IDS,GM,Minf,WNP1,P)
 Implicit None
!*******************************************************************************************
 Intent(In)::Dim,NFW1,NFW2,NC,NP,NF,X,Y,IDS,GM,Minf,WNP1,P

 Integer::Dim,I,J,P1,P2,P3,P4,ME,NP,NC,NFW1,NFW2,NF
 Real(8)::Minf,GM,CP,CCP,Xm
 Integer,Dimension(1:4,1:Dim)::Corn,IDS
 Real(8),Dimension(1:Dim)::X,Y,P
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*******************************************************************************************	
!Part 1:
 Open(1,File='Contours.Plt')
 Open(2,File='CP.Plt')
 Open(3,File='SolutionData.txt')

!Part 1:
 Call EbasedToCbased(Dim,NC,NF,IDS,Corn)


!Part 1:
 Write(1,*) 'TITLE="All Data"'
 Write(1,*) 'Variables="X","Y","Ro","U","V","P"'
 Write(1,*) 'ZONE N=',NP,' E=',NC,' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK VARLOCATION=([3-6]=CELLCENTERED)'

 Do J=1,NP
	Write(1,*) X(J)
 End Do
 Do J=1,NP
	Write(1,*) Y(J) 
 End Do

 Do J=1,NC
	Write(1,*) WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) WNP1(2,J)/WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) WNP1(3,J)/WNP1(1,J)
 End Do
 Do J=1,NC
	Write(1,*) P(J)
 End Do
 
 Do I=1,NC
    P1 = Corn(1,I) 
    P2 = Corn(2,I) 
    P3 = Corn(3,I)
    P4 = Corn(4,I)
    if(P4==0) P4=P3

	Write(1,*) P1,P2,P3,P4
 End Do



!Part 1:
 CCP = 0.5*GM*Minf*Minf
 WRITE(2,*)'VARIABLES="X","CP"'  
 WRITE(2,*)'ZONE'

 DO J=NFW1+1,NFW2
    ME  = IDS(1,J)
    P1  = IDS(3,J)
    P2  = IDS(4,J)

    Xm = 0.5*( X(P1)+X(P2) )
    CP = (GM*P(ME)-1.0)/CCP

	Write(2,*) Xm,-CP
 end do



!Part 1:
 Do I=1,NC
	Write(3,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I)
 End Do


 Close(1)
 Close(2)
 Close(3)
!*******************************************************************************************
 End
!###########################################################################################

!*******************************************************************************************
 Subroutine EbasedToCbased(Dim,NC,NF,IDS,Corn)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF,IDS 
 Intent(Out  )::Corn

 Integer::Dim,J,I,NC,NF,ME,NE,P1,P2 
 Integer,Dimension(1:4,1:Dim)::IDS,Corn
 Integer,Dimension(1:Dim)::NCorn
!*******************************************************************************************

!Part 6:
 Do J=1,NC
	Corn(1,J)=0
	Corn(2,J)=0
	Corn(3,J)=0
	Corn(4,J)=0
    NCorn(J)=0
 End Do

!Part 6:
 Do J=1,NF
    ME = IDS(1,J)
    NE = IDS(2,J)
    P1 = IDS(3,J)
	P2 = IDS(4,J)

	Do I=1,NCorn(ME)
	   IF( P1==Corn(I,ME) ) goto 1
    End Do
    NCorn(ME) = NCorn(ME) + 1
	Corn( NCorn(ME),ME ) = P1

1	Do I=1,NCorn(ME)
	   IF( P2==Corn(I,ME) ) goto 2
    End Do
    NCorn(ME) = NCorn(ME) + 1
	Corn( NCorn(ME), ME) = P2

2   IF(NE==0)goto 4

	Do I=1,NCorn(NE)
	   IF( P1==Corn(I,NE) ) goto 3
    End Do
    NCorn(NE) = NCorn(NE) + 1
	Corn(NCorn(NE), NE ) = P1

3	Do I=1,NCorn(NE)
	   IF( P2==Corn(I,NE) ) goto 4
    End Do
    NCorn(NE) = NCorn(NE) + 1
	Corn(NCorn(NE), NE ) = P2  
 
4 End Do
!*******************************************************************************************
 End
!###########################################################################################
 Subroutine Read_2DMesh1(Dim,NP1,NC1,NF_1,NR1,NFR1,BC1,IDS1,X1,Y1,NOL,FTC)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NOL,FTC
 Intent(Out  )::NP1,NC1,NF_1,NR1,NFR1,BC1,IDS1,X1,Y1

 Integer::Dim,I,J,J1,JJ,NP1,NC1,NF_1,NR1,SFace,FaceType,MeshDim,NOL,FTC
 Integer,Dimension(1:100)::NFR1,BC1
 Integer,Dimension(1:4,1:Dim)::IDS1
 Real(8),Dimension(1:Dim)::X1,Y1
!*******************************************************************************************
!Part 1:
 IF(NOL==1 .AND. FTC==1)then
   Open(1,File='Mesh1.Txt')
 else IF(NOL==1 .AND. FTC==0)then
   Open(1,File='Mesh2.Txt')
 else if(NOL==2 .AND. FTC==1)then
   Open(1,File='Mesh2.Txt')
 else if(NOL==2 .AND. FTC==0)then
   Open(1,File='Mesh3.Txt')
 end if

!Part 2:
 Read(1,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(1,*) NP1    

!Part 4:
 Read(1,*) NC1

!Part 5:
 Read(1,*) NF_1

!Part 6:
 Read(1,*) NR1

!Part 7:
 Read(1,*)
 Do J=1,NR1
    Read(1,*) NFR1(J) , BC1(J)
 End Do 

!Part 8:
 Read(1,*)
 Do J=1,NF_1  
    Read(1,*) FaceType,IDS1(1,J),IDS1(2,J),IDS1(3,J),IDS1(4,J)
 End Do

!Part 9:
 Read(1,*)        
 Do J=1,NP1
    Read(1,*) X1(J),Y1(J)
 End Do

 Close(1)
!*******************************************************************************************
End
!###########################################################################################
Subroutine Read_2DMesh2(Dim,NP2,NC2,NF_2,NR2,NFR2,BC2,IDS2,X2,Y2,NOL,FTC)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NOL,FTC
 Intent(Out  )::NP2,NC2,NF_2,NR2,NFR2,BC2,IDS2,X2,Y2

 Integer::Dim,I,J,J1,JJ,NP2,NC2,NF_2,NR2,SFace,FaceType,MeshDim,NOL,FTC
 Integer,Dimension(1:100)::NFR2,BC2
 Integer,Dimension(1:4,1:Dim)::IDS2
 Real(8),Dimension(1:Dim)::X2,Y2
!*******************************************************************************************
!Part 1:
 IF(NOL==1 .AND. FTC==1)then
   Open(2,File='Mesh2.Txt')
 else IF(NOL==1 .AND. FTC==0)then
   Open(2,File='Mesh1.Txt')
 else if(NOL==2 .AND. FTC==1)then
   Open(2,File='Mesh3.Txt')
 else if(NOL==2 .AND. FTC==0)then
   Open(2,File='Mesh2.Txt')
 end if

!Part 2:
 Read(2,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(2,*) NP2    

!Part 4:
 Read(2,*) NC2

!Part 5:
 Read(2,*) NF_2

!Part 6:
 Read(2,*) NR2

!Part 7:
 Read(2,*)
 Do J=1,NR2
    Read(2,*) NFR2(J) , BC2(J)
 End Do 

!Part 8:
 Read(2,*)
 Do J=1,NF_2  
    Read(2,*) FaceType,IDS2(1,J),IDS2(2,J),IDS2(3,J),IDS2(4,J)
 End Do

!Part 9:
 Read(2,*)        
 Do J=1,NP2
    Read(2,*) X2(J),Y2(J)
 End Do

 Close(2)
!*******************************************************************************************
End
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Find Index of Points Constructing the cell                           //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenMesh@chmail.ir                            //!
!// Doc ID: MC5F088F1                                                                    //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Edge_To_Cell(Dim,NF,NC,IDS,Corn)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS
 Intent(Inout)::Corn
 
 Integer::Dim,NF,NC,ME,NE,I,J,J1,J2,E,E1,E2,E3,P1_E2,P2_E1,P
 Integer,Dimension(1:4, 1:Dim)::IDS
 Integer,Dimension(1:4, 1:Dim)::CELL_EDGE,Corn
 Integer,Dimension(1:Dim)::NCELL_EDGE
!*******************************************************************************************
!Part 1:
 Do I=1,NC
    NCELL_EDGE(I) = 0   
	Corn(1,I)     = 0  
	Corn(2,I)     = 0  
	Corn(3,I)     = 0  
	Corn(4,I)     = 0
 End Do
 
!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1, I)
    NE = IDS(2, I)
	
   !Part 4:
	NCELL_EDGE(ME) = NCELL_EDGE(ME) + 1
    CELL_EDGE(NCELL_EDGE(ME),ME)=I

   !Part 5:
    IF(NE/=0)Then
	 NCELL_EDGE(NE) = NCELL_EDGE(NE) + 1
     CELL_EDGE(NCELL_EDGE(NE),NE)=-I
    EndIF

 End Do

!Part 6:
 Do I=1,NC

   !Part 7:
    Do J1=1,NCELL_EDGE(I)
       E1 = CELL_EDGE(J1,I)

      !Part 8:
       IF( E1>0 )Then
        P2_E1 = IDS(4,E1)
       Else
        P2_E1 = IDS(3,-E1)
	   EndIF

      !Part 9:
	   Do J2=J1+1,NCELL_EDGE(I)
          E2 = CELL_EDGE(J2,I)

         !Part 10:
          IF( E2>0 )Then
           P1_E2 = IDS(3,E2)
          Else
           P1_E2 = IDS(4,-E2)
          EndIF

         !Part 11:
          IF( P2_E1==P1_E2 )Then
	       E                 = CELL_EDGE(J1+1,I)
	       CELL_EDGE(J1+1,I) = CELL_EDGE(J2,I)
           CELL_EDGE(J2,I)   = E
          EndIF

       End Do

    End Do

 End Do 

!Part 12:
 Do I=1,NC
    Do J=1,NCELL_EDGE(I)

       E = CELL_EDGE(J,I)

       IF( E>0 )Then
        P = IDS(3,E)
       Else
        P = IDS(4,-E)
       EndIF

       Corn(J,I) = P

    End Do
 End Do
!*******************************************************************************************
 End
!###########################################################################################
 !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Calculate Residual of Mass Equation                                   //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F004F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine ResMass1(Dim,NC,WNP1,WN,DT,Rm)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,WNP1,WN ,DT
 Intent(Out  )::Rm

 Integer::Dim,I,NC
 Real(8)::Rm
 Real(8),Dimension(1:4,1:Dim)::Wnp1,WN
 Real(8),Dimension(1:Dim)::DT
!*******************************************************************************************	
!Part 1:
 Rm = 0.0

!Part 2:
 DO I=1,NC
    Rm = Rm + Abs( Wnp1(1,I)-WN(1,I) ) / WN(1,I)
 End Do

!Part 3:
 Rm = Rm / NC

!Part 4:
 Rm = Dlog10(Rm+1.0E-16)
!*******************************************************************************************
 End
!###########################################################################################
 Subroutine Mesh_Info(Dim,IDS,X,Y,CornC,NeibC,NC,NF)
Implicit None
!*******************************************************************************************
Intent(In   )::Dim,NC,NF,IDS
Intent(Out  )::CornC,NeibC

Integer::Dim  
Integer::I,J,ME,NE,P1,P2,counter,NC,NF
Integer,Dimension(1:3,1:Dim)::CornC,NeibC
Real(8),Dimension(1:Dim)::X,Y
Integer,Dimension(1:4,1:Dim)::IDS
!*******************************************************************************************    
Do I=1,NC
    counter=0
    Do J=1,NF
        ME=IDS(1,J)
        NE=IDS(2,J)
        P1=IDS(3,J)
        P2=IDS(4,J)
        IF (I==ME .AND. counter==0) Then
            CornC(1,I)=P1
            CornC(2,I)=P2
            NeibC(3,I)=NE
            counter=counter+1
            cycle
        End IF
            
        IF (I==NE .AND. counter==0) Then
            CornC(1,I)=P1
            CornC(2,I)=P2
            NeibC(3,I)=ME
            counter=counter+1
            cycle
        End IF
            
        IF (I==ME) Then
                
            IF((CornC(3,I)==P1 .AND. CornC(2,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(2,I)==P1)) Then
                counter=counter+1
                NeibC(1,I)=NE
                cycle
            End IF
                
            IF((CornC(3,I)==P1 .AND. CornC(1,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(1,I)==P1)) Then
                counter=counter+1
                NeibC(2,I)=NE
                cycle
            End IF
                
            IF (CornC(1,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(2,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(1,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(2,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(1,I)=NE;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(1,I)=NE;
                Cycle;
            End IF
                
        End IF
            
            
        IF (I==NE) Then
                
            IF((CornC(3,I)==P1 .AND. CornC(2,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(2,I)==P1)) Then
                counter=counter+1
                NeibC(1,I)=ME
                cycle
            End IF
                
            IF((CornC(3,I)==P1 .AND. CornC(1,I)==P2) .OR. (CornC(3,I)==P2 .AND. CornC(1,I)==P1)) Then
                counter=counter+1
                NeibC(2,I)=ME
                cycle
            End IF
                
            IF (CornC(1,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(2,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(1,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(2,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P1) Then
                counter=counter+1;
                CornC(3,I)=P2;
                NeibC(1,I)=ME;
                Cycle;
            End IF
                
            IF (CornC(2,I)==P2) Then
                counter=counter+1;
                CornC(3,I)=P1;
                NeibC(1,I)=ME;
                Cycle;
            End IF
        End IF
                    
    End Do
End Do
           
!Do I=1,NC
!    print*,CornC(1,I),CornC(2,I),CornC(3,I),NeibC(1,I),NeibC(2,I),NeibC(3,I)
!End Do

End


