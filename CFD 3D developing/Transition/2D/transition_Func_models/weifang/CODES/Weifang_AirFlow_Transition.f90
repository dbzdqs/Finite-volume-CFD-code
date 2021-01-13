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
 Integer,Parameter::Dim=230000
!===============================

 Integer::I,J,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NRKS,NWrite,Init,Ncyc,NS , NN
 Real(8)::GM,Co,ALF,R0,P0,T0,C0,U0,V0,E0,B0,ERmx,CFLx,Minf,Rinf,Tt,MR,PrL,PrT,Rm,U,V,RKco,Mu0,&
          Mut0,Temp,TUinf,Mutinf,e,CD,CL,ALFA_GRAPH,ALF_MAX,ALF_STEP,INV_INIT
 Integer,Dimension(1:100)::NFR,BC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1,WC,Con,Dif
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,A,NX,NY,DA,DW,DT,P,Mu,Mut,DUX,DUY,DVX,DVY,DTX,DTY,Taukk
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Wnt,Wbt
 Real(8),Dimension(1:4,1:Dim)::Limit
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
 Real(8),Dimension(1:3,1:Dim)::R_LSQ
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_Setting_V3(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,TUinf,Mutinf,e,ALFA_GRAPH,ALF_MAX,ALF_STEP,NN,INV_INIT)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 4:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 
!Part 6:
 Call InitMeanFlow(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,T0,Mu0,B0,WNP1)
 
 
 !*************************
 
 IF (INV_INIT==0.0) goto 198
 
 !*************************
 
   
!Part 7: 
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

!Part 8:
 
 Call BC_InvisWall(Dim,NFW1,NFW2,NX,NY,DA,IDS,GM,WNP1,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow_V1(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 9:
  Call FirstOrd_Gradient(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1)
 !Call LSQ_GEO_COEFF(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,XC,YC,X,Y,R_LSQ)
 !Call LSQ_GRADIENT2(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y,R_LSQ)
 ! Call LSQ_GRADIENT(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y)

!Part 10:
 Call LIMITER_v2(Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,Limit,X,Y,NF,e)

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
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)


   !Part 16:
    Do NS=1,2
   
      !Part 17:
	 RKco=1.0/(2-NS+1)
       
     !  Do NS=1,1
   
      !Part 17:
	 !  RKco=1.0/(1-NS+1)

      !Part 18:
	   
    !  Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
      Call ConMeanFlow_AUSM_PlusUP_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,GWNP1,Xc,Yc,Limit,Con,X,Y)
    
   
      !Part 19:
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

      !Part 20:
       
       Call BC_InvisWall(Dim,NFW1,NFW2,NX,NY,DA,IDS,GM,WNP1,P,WB)
       !Call BC_Wall_Inv_Press2nd(Dim,NFW1,NFW2,NX,NY,DA,IDS,GM,WNP1,P,WB,X,Y,XC,YC,GWNP1,Limit)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow_V1(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

      !Part 21: 
       Call FirstOrd_Gradient(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1)
      !Call LSQ_GEO_COEFF(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,XC,YC,X,Y,R_LSQ)
      !Call LSQ_GRADIENT2(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y,R_LSQ)
    ! Call LSQ_GRADIENT(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y)
      !Part 22: 
       Call LIMITER_v2(Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,Limit,X,Y,NF,e)

     End Do !Ns	
     
	!Part 23:
 	Call ResMass(Dim,NC,WNP1,WC,DT,Rm)
     Print*,Ncyc,Rm  ,'INVISCID RUN'
     
     IF( Mod(Ncyc,2000)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
	  
      Open(30,File='SolutionData.txt')
      Do I=1,NC
	     Write(30,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I)
      End Do
      close(30)
      
    end if

 End Do !Do While

Print*,'**********************END OF INITIALIZATION*************************'


!***************************************************************************************

198 Do J=1,NC
     
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
 Call BC_InFlow_V1(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

 
 Call FirstOrd_Gradient(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1)

! Call LSQ_GEO_COEFF(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,XC,YC,X,Y,R_LSQ)
 !Call LSQ_GRADIENT2(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1,XC,YC,X,Y,R_LSQ)
 
 
 Call Limiter_V2(Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,Limit,X,Y,NF,e)
 
 
!Part 9:
 Call Gamma_TraSST_Init(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,Dw,INW,R0,Minf,Rinf,Wnt,Wntp1,Mut,TUinf,Mutinf,Init,WNP1)  



!mut=0.1

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
	   RKco=1.0/(NRKS-NS+1)

      !Part 18:
	   !Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
        Call ConMeanFlow_AUSM_PlusUP_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,GWNP1,Xc,Yc,Limit,Con,X,Y)

      !Part 19:
       Call VelTemp_GradFace(Dim,NC,NF1,NF2,NFW1,NFW2,NF,NP,IDS,X,Y,Xc,Yc,WNP1,WB,GM,P,DUX,DUY,DVX,DVY,DTX,DTY)

      !Part 20:
	   Call DifMeanFlow_Turb_Rk(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,GM,PrL,PrT,NX,NY,MR,Mu,Mut,WNP1,Wntp1,WB,&
	                            DUX,DUY,DVX,DVY,DTX,DTY,Dif)

           
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
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))

         !Part 23:
          Temp = GM*P(J)/WNP1(1,J) 
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         

       End Do

      !Part 24: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       !Call BC_Wall_Visc_Press2nd(Dim,NFW1,NFW2,IDS,GM,P,WB,WNP1,X,Y,XC,YC,GWNP1,Limit)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow_V1(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
       
       Call FirstOrd_Gradient(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,WNP1,WB,P,GWNP1)

      !Part 22: 
       Call Limiter_V2(Dim,NC,NF1,NF2,IDS,GM,XC,YC,WNP1,P,WB,GWNP1,Limit,X,Y,NF,e)

     End Do !Ns	

    !Part 25:
        
     Call Weifang_Transition_SST(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,Minf,Rinf,MR,R0,&
            NRKS,Mu0,Wb,Wbt,Wnp1,Wnt,Wntp1,Mu,Mut, NN,TUinf,Mutinf)

     !===========================================================================/
     Call PressLiftDragCo(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,CL,CD)
     Open(1110,File='CD_Iter.Plt')
 	 Open(1120,File='CL_Iter.Plt')
     Write(1110,*) CD
     Write(1120,*) CL
    !============================================================================/
    !Part 26:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Rm)
     Print*,Ncyc,Rm,'TRANSITION_TURBULENT_RUN'
     
     
      IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
	  Call Write_ResultsV3(Dim,NFW1,NFW2,NF,NC,NP,IDS,X,Y,Xc,Yc,GM,Minf,Rinf,WNP1,P,Mu,Mut,DUY,ALF,Wntp1)
      !Call Write_Results_V2(Dim,NFW1,NFW2,NF,NC,NP,IDS,X,Y,Xc,Yc,GM,Minf,Rinf,WNP1,P,Mu,Mut,DUY,ALF)

	 End If

 End Do !Do While
!============================================================================== \
 
 
 If (ALFA_GRAPH==1.0) then
 
 NN=NN+1
 Call Write_ResultsV4(Dim,NFW1,NFW2,NF,NC,NP,IDS,X,Y,Xc,Yc,GM,Minf,Rinf,WNP1,P,Mu,Mut,DUY,ALF,Wntp1,NN)
   
     Open(1200,File='CD_ALFA.Plt')
 	 Open(1210,File='CL_ALFA.Plt')
     Write(1200,*) ALF,CD
     Write(1210,*) ALF,CL
     
     ALF=ALF+ALF_STEP
     
   IF(ALF<=ALF_MAX)  goto 198
     
 end if
     
  
 !==============================================================================\
  Call Write_ResultsV3(Dim,NFW1,NFW2,NF,NC,NP,IDS,X,Y,Xc,Yc,GM,Minf,Rinf,WNP1,P,Mu,Mut,DUY,ALF,Wntp1)
  !Call Write_Results_V2(Dim,NFW1,NFW2,NF,NC,NP,IDS,X,Y,Xc,Yc,GM,Minf,Rinf,WNP1,P,Mu,Mut,DUY,ALF)
 !=============================================================================/
   Open(1130,File='CD-CL.txt')
   Write(1130,*) 'CD equales to:' , CD
   Write(1130,*) 'CL equales to:' , CL
   close(1130)
 !=============================================================================/
   
 PRINT*,'*********END OF CALCULATION**********'
 PAUSE
!*******************************************************************************************
 End 
!###########################################################################################


