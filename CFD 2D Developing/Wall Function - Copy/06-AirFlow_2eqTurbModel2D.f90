!//////////////////////////////////////////////////////////////////////////////////////////!
!// AirFlow_Turb                                                                         //!
!// Date :         Febreury/2/2015                                                       //!
!// Developed by : M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                          //!
!//                                                                                      //!                                                                      //
!// A Turbulent 2D Flow Solver                                                           //!
!// Features: 1- 2D                                                                      //!
!//           2- Using Unstructured Mesh                                                 //!
!//           3- Edge Based Data Structured                                              //!
!//           4- Cell Center Conrol Volume                                               //!
!//           5- Laminar Flow                                                            //!
!//           6- Convection Terms is Discritized by AUSM Scheme                          //!
!//           7- Transient Term is Discritized by Runge-Kutta Explicit                   //!
!//           8- Prandtl Tulbulence Model                                                //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Program AirFlow_Turb
 Implicit None
!===============================
 Integer,Parameter::Dim=200000
!===============================

 Integer::I,J,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NRKS,NWrite,Init,Ncyc,NS  
 Real(8)::GM,Co,ALF,R0,P0,T0,C0,U0,V0,E0,B0,ERmx,CFLx,Minf,Rinf,Tt,MR,PrL,PrT,Rm,U,V,RKco,Mu0,&
          Mut0,Temp
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,Dif
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,A,NX,NY,DA,DW,DT,P,Mu,Mut,DUX,DUY,DVX,DVY,DTX,DTY,Taukk
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:Dim)::INW
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell
 Integer,Dimension(1:Dim)::NEdgOfCell
 Integer,Dimension(1:4, 1:Dim)::Corn
 Real(8)::K2,K4,TotTime,Time_Coe
 Real(8),Dimension(1:5)::RKJ
 Real(8),Dimension(1:3,1:Dim)::iwf
 Integer,Dimension(1:2, 1:Dim)::nn
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
 Call PointOfCell(Dim,NC,InxEdgOfCell,NEdgOfCell,IDS,Corn)
 
 Call WallDist(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,INW,DW)
 
!Part 4:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 5:
 Call InitMeanFlow(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,T0,Mu0,B0,WNP1)

 Do J=1,NC
     
   !Part 6: 
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

   !Part 7:
    Temp = GM*P(J)/WNP1(1,J) 
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         

 End Do

!Part 8:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 9:for all ke turb model
 !Call Ke_Init(Dim,NC,MR,WTNP1,Mut)
 Call Ke_Init2(Dim,NC,MR,iwf,WTNP1,Mut,Mu,WNP1) 
 !!!Call Kw_Init(Dim,NC,MR,WTNP1,Mut) 
 !call Ke_Init2(Dim,NC,NF,NFW1,NFW2,IDS,X,Y,Xc,Yc,MR,DW,INW,IWF,WTNP1,Mut,NN,wnp1,mu) 
      
!Part 10:
 
!Part 11:
 Ncyc = 0
 Rm   = 10.0

!Part 12:
 Do While(Rm > ERmx)

   !Part 13:
    Ncyc=Ncyc+1
     
   !Part 14:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   !Part 15:
	Call TimSTP_Turb(Dim,NC,NF,NF1,NF2,IDS,NX,NY,DA,A,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)

   !Part 16:
    Do NS=1,NRKS
   
      !Part 17:
	   RKco=RKJ(NS)

      !Part 18:
	   Call ConMeanFlow_ScalarDiss(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,K2,K4,WNP1,WB,P,Con)

      !Part 19:
       Call VelTemp_GradFace(Dim,NC,NF1,NF2,NFW1,NFW2,NF,NP,IDS,X,Y,Xc,Yc,WNP1,WB,GM,P,DUX,DUY,DVX,DVY,DTX,DTY)

      !Part 20:
       call DifMeanFlow_Turb_wf(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,GM,PrL,PrT,NX,NY,MR,Mu,Mut,&
                                WNP1,WTNP1,WB,DUX,DUY,DVX,DVY,DTX,DTY,Dif,iwf,DW)
       
      !Part 21:
       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co*( Con(1,J)          )
		  WNP1(2,J) = WC(2,J) - Co*( Con(2,J)+Dif(2,J) )
          WNP1(3,J) = WC(3,J) - Co*( Con(3,J)+Dif(3,J) )
          WNP1(4,J) = WC(4,J) - Co*( Con(4,J)+Dif(4,J) )

         !Part 22:
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V)-0.5*WTNP1(1,J)) !

         !Part 23:
          Temp = GM*P(J)/WNP1(1,J) 
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         

       End Do

      !Part 24: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

    End Do !Ns	

    !Part 25:
     !!!Call KeChien_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
     !!!                  NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,WTNP1,Mut)

     !!!Call KeLB_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
     !!!                 NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,P,GM,DUY,WTNP1,Mut)
     
     !!!Call KeYang_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
     !!!                    NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,WTNP1,Mut)
     
     !!!Call KwWilcox_Main(Dim,NC,NP,DUY,INW,IDS,MR,NRKS,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,DA,&
     !!!                   NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,NX,NY,XC,YC,DW,A,Mu,DT,Mut)
     
     !!!Call KwBredberg_Main(Dim,NC,NP,DUY,INW,IDS,MR,NRKS,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,&
     !!!                       NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,NX,NY,XC,&
					!!!		YC,DW,A,Mu,DT,Mut)
     
     !!!Call KWSST_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
     !!!                NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
     
     !!!Call KWSST_Sust_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
     !!!                NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
    
     !!!Call KWSST_V_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
     !!!                NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,Mu0,Wb,WNP1,Mu,WTNP1,Mut)    
    
     call Ke_Mainw(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                         NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,WTNP1,Mut,iwf,RKJ,p)

    !Part 26:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,Ncyc,Rm

    !Part 27:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
      Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
      Call Write_CF(Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu)
      Call Write_ConservativeVariables(Dim,NC,WNP1)
	  Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mu,"MolVis.plt")
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mut,"TubVis.plt")
	 End If

 End Do !Do While
!*******************************************************************************************
 End 
!###########################################################################################


