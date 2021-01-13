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
 Program Multigrid
 Implicit None
!===============================
 Integer,Parameter::Dim=120000
!===============================
 Integer::I,J,NS
 Real(8)::U,V,RKco,Co
 
 Integer::NC
 Integer::NP
 Integer::NF
 Integer::NF1,NF2
 Integer::NFW1,NFW2
 Integer::NFF1,NFF2
 Integer::NFI1,NFI2
 Integer::NFS1,NFS2
 Integer::NFO1,NFO2
 Integer::NFIF1,NFIF2
 Integer::NRKS
 Integer::NWrite
 Integer::NR
 Integer::Init
 Integer::Ncyc
 Real(8)::GM
 Real(8)::ALF
 Real(8)::R0
 Real(8)::P0
 Real(8)::C0
 Real(8)::U0
 Real(8)::V0
 Real(8)::E0
 Real(8)::H0
 Real(8)::ERmx
 Real(8)::CFLx
 Real(8)::Minf
 Real(8)::Rm
 Integer,Dimension(1:100)::NFR
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:4,1:Dim)::WC
 Real(8),Dimension(1:4,1:Dim)::Con
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:Dim)::XC,YC
 Real(8),Dimension(1:Dim)::A
 Real(8),Dimension(1:Dim)::NX,NY
 Real(8),Dimension(1:Dim)::DA
 Real(8),Dimension(1:Dim)::DT
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell
 Integer,Dimension(1:Dim)::NEdgOfCell
 Integer,Dimension(1:4, 1:Dim)::Corn
 Real(8)::Rinf
 Real(8)::Tt
 Real(8)::K2,K4
 Real(8)::TotTime
 Real(8)::Time_Coe
 Real(8),Dimension(1:5)::RKJ
 
 Integer::counter,C,NSN,NMF,NE,ME,NOL,FTC,NON
 Real(8)::DTmin,ERmx1,epsilon,ER1,ER2,ER3,ER4
 Integer,Dimension(1:3,1:Dim)::Corn1,Corn2,CornC,NeibC
 Real(8),Dimension(1:4,1:Dim)::RES,RESM0,RESM1,RESM2,Error,SolM1,WC1,FFM2
 Real(8),Dimension(1:Dim)::X2,Y2,A1,X1,Y1
 real::start,finish
!************************************************** Main **************************************************
 call cpu_time(start)

!Part 1:
 NMF=1
 Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF) 

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
 Call PointOfCell(Dim,NC,InxEdgOfCell,NEdgOfCell,IDS,Corn)
 
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
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
 
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
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

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
!!!Call Interpolation(Dim,WNP1,RES,Error,NOL,FTC)

 !Part 20:
 NMF=2
 Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF) 
 
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
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
 
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
ERmx1=-5.5

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
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

    End Do !Ns	

    !Part 39:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
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
 !!Call Interpolation(Dim,WNP1,RES,Error,NOL,FTC)
 
 !Part 42:
 NMF=1
 Call Read_2DMeshMG(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,NMF) 
 
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
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

 
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
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

    End Do !Ns	
    
    Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,Ncyc,Rm

     !Part 49:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
        Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
        Call Write_ConservativeVariables(Dim,NC,WNP1)
        Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
	 End If
    
 End Do !while
 
Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
Call Write_ConservativeVariables(Dim,NC,WNP1)
Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
 call cpu_time(finish)
 print *, finish-start
    
 pause
!**********************************************************************************************************
 End 
!##########################################################################################################
