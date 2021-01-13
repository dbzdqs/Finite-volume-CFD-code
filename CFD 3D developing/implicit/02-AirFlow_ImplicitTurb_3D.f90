!//////////////////////////////////////////////////////////////////////////////////////////!
!// AirFlow_Turb3D                                                                       //!
!// Date :         Febreury/2/2015                                                       //!
!// Developed by : B. Jodeiri, Iran, Tehran, b.jodeiri@ut.ac.ir                          //!
!//                                                                                      //!
!// An Inviscid 3D Flow Software                                                         //!
!// Features: 1- 3D                                                                      //!
!//           2- Using Unstructured Mesh                                                 //!
!//           3- Edge Based Data Structured                                              //!
!//           4- Cell Based Conrol Volume                                                //!
!//           5- Turb Flow                                                               //!
!//           6- Convection Terms is Discritized by AUSM Scheme                          //!
!//           7- Transient Term is Discritized by Runge-Kutta Explicit                   //!
!//                                                                                      //!
!// The Program Is Available Through The Website: WWW.IRCSE.IR                           //!
!// It May Be Copied, Modified, And Redistributed For Non-Commercial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Program AirFlow_Turb3D
 Implicit None
!===============================
 Integer,Parameter::Dim=500000
!===============================

 Integer::I,J,NC,NP,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
          NR,NRKS,NWrite,Init,Ncyc,NS ,j1,j2,j3 ,j4,j5,cell,Naverage,NFEmp
 Real(8)::GM,Co,ALF,R0,P0,C0,U0,V0,W0,E0,ERmx,CFLx,Minf,Tt,B0,Rm,U,V,W,RKco,Temp,Dmin,Rinf,MR,&
          PrL,PrT,T0,Mu0,Time,Mut0,DT_min,ERmxStdy,TotTime,ERmxUnStdy,Time_Coe,K2,K4
 Integer,Dimension(1:100)::NFR,BC
 
 Integer,allocatable,Dimension(:,:)::IDS
 Real(8),allocatable,Dimension(:,:)::WNP1
 Real(8),allocatable,Dimension(:,:)::WC
 Real(8),allocatable,Dimension(:,:)::Con
 Real(8),allocatable,Dimension(:,:)::Dif
 Real(8),allocatable,Dimension(:)::X,Y,Z
 Real(8),allocatable,Dimension(:)::XC,YC,ZC
 Real(8),allocatable,Dimension(:)::Vol
 Real(8),allocatable,Dimension(:)::NX,NY,NZ
 Real(8),allocatable,Dimension(:)::DA
 Real(8),allocatable,Dimension(:)::DT
 Real(8),allocatable,Dimension(:)::P
 Real(8),allocatable,Dimension(:)::Mu
 Real(8),allocatable,Dimension(:)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX
 Real(8),allocatable,Dimension(:)::DTY,DTZ
 Real(8),allocatable,Dimension(:)::TurbQ
 Real(8),allocatable,Dimension(:,:)::WB
 Integer,allocatable,Dimension(:)::FaceType
 Integer,allocatable,Dimension(:)::N_CORN
 Integer,allocatable,Dimension(:,:)::IFace_Cell
 Integer,allocatable,Dimension(:)::NFace_Cell
 Integer,allocatable,Dimension(:,:)::Corn
 Integer,allocatable,Dimension(:)::INW
 Real(8),allocatable,Dimension(:)::DW
 Real(8),allocatable,Dimension(:)::Mut
 Real(8),allocatable,Dimension(:,:)::WTNP1
 Real(8),allocatable,Dimension(:,:)::WTB
 Real(8),allocatable,Dimension(:)::SumCF
 Real(8),allocatable,Dimension(:)::SumCP
 Real(8),allocatable,Dimension(:)::sum_umean
 Real(8),allocatable,Dimension(:)::sum_vmean
 Real(8),allocatable,Dimension(:)::sum_uvmean
 Real(8),allocatable,Dimension(:)::sum_u2mean
 Real(8),allocatable,Dimension(:)::sum_v2mean
 
 Real(8),Dimension(1:5,1:Dim)::Limit
 Real(8),Dimension(1:3,1:5,1:Dim)::GWNP1
 
 Real(8),Dimension(1:5)::RKJ
 Real(8)::Kinf
 Real(8)::Omegainf
 Real(8)::Gamainf
 Real(8)::Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0
 Real(8)::TUinf
 
 

 Real(8)::sor
 
 !!!Real(8),Dimension(1:5,1:Dim)::DQ
 !!!Real(8),Dimension(1:5,1:Dim)::Res
 !!!Integer,Dimension(1:6, 1:Dim)::INeib
 !!!Integer,Dimension(1:10,1:Dim)::InonzeroCell
 !!!Integer,Dimension(1:10,1:Dim)::InonzeroFace
 !!!Integer,Dimension(1:Dim)::NnonzeroCell
 !!!Integer,Dimension(1:Dim)::NnonzeroFace
 !!!Integer,Dimension(1:Dim)::NUM_NONZERO
 
 Real(8),allocatable,Dimension(:,:)::DQ
 Real(8),allocatable,Dimension(:,:)::Res
 Integer,allocatable,Dimension(:,:)::INeib
 Integer,allocatable,Dimension(:,:)::InonzeroCell
 Integer,allocatable,Dimension(:,:)::InonzeroFace
 Integer,allocatable,Dimension(:)::NnonzeroCell
 Integer,allocatable,Dimension(:)::NnonzeroFace
 Integer,allocatable,Dimension(:)::NUM_NONZERO
!***************************************** Main ********************************************
!Part 1:
 allocate( IDS(1:6,1:Dim)        )
 allocate( WNP1(1:5,1:Dim)       )
 allocate( WC(1:5,1:Dim)         )
 allocate( Con(1:5,1:Dim)        )
 allocate( X(1:Dim)              )
 allocate( Y(1:Dim)              )
 allocate( Z(1:Dim)              )
 allocate( XC(1:Dim)             )
 allocate( YC(1:Dim)             )
 allocate( ZC(1:Dim)             )
 allocate( Vol(1:Dim)            )
 allocate( NX(1:Dim)             )
 allocate( NY(1:Dim)             )
 allocate( NZ(1:Dim)             )
 allocate( DA(1:Dim)             )
 allocate( TurbQ(1:Dim)          )
 allocate( WB(1:6,1:Dim)         )
 allocate( FaceType(1:Dim)       )
 allocate( N_CORN(1:Dim)         )
 allocate( IFace_Cell(1:6,1:Dim) )
 allocate( NFace_Cell(1:Dim)     )
 allocate( Corn(1:8,1:Dim)       )
 allocate( DT(1:Dim)             )
 allocate( P(1:Dim)              )
 
!Part 2:
 allocate( INW(1:Dim)       )
 allocate( DW(1:Dim)        ) 
 allocate( Dif(1:5,1:Dim)   )
 allocate( Mu(1:Dim)        )
 allocate( Mut(1:Dim)       )
 allocate( DUX(1:Dim)       )
 allocate( DUY(1:Dim)       )
 allocate( DUZ(1:Dim)       )
 allocate( DVX(1:Dim)       )
 allocate( DVY(1:Dim)       )
 allocate( DVZ(1:Dim)       )
 allocate( DWX(1:Dim)       )
 allocate( DWY(1:Dim)       )
 allocate( DWZ(1:Dim)       )
 allocate( DTX(1:Dim)       )
 allocate( DTY(1:Dim)       )
 allocate( DTZ(1:Dim)       )
 allocate( WTNP1(1:2,1:Dim) )
 allocate( WTB(1:2,1:Dim)   )
 
 allocate( DQ(1:5,1:Dim)            )
 allocate( Res(1:5,1:Dim)           )
 allocate( INeib(1:6, 1:Dim)        )
 allocate( InonzeroCell(1:10,1:Dim) )
 allocate( InonzeroFace(1:10,1:Dim) )
 allocate( NnonzeroCell(1:Dim)      )
 allocate( NnonzeroFace(1:Dim)      )
 allocate( NUM_NONZERO(1:Dim)       )
 
!Part 3:
 allocate( SumCF(1:Dim)      )
 allocate( SumCP(1:Dim)      )
 allocate( sum_umean(1:Dim)  )
 allocate( sum_vmean(1:Dim)  )
 allocate( sum_uvmean(1:Dim) )
 allocate( sum_u2mean(1:Dim) )
 allocate( sum_v2mean(1:Dim) )
 SumCF=0.0;SumCP=0.0;sum_umean=0.0;sum_vmean=0.0;sum_uvmean=0.0;sum_u2mean=0.0;sum_v2mean=0.0
 
 
 
!Part 4:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)
 
!Part 5:
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"Imesh.plt")
 print*,'read mesh finished'
 
!Part 6:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)

!Part 7:
 Call MeshBC3D(Dim,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 8:
 Call FaceOfCell(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 
 Call NeibOfCell3D(Dim,NC,IDS,NFace_Cell,IFace_Cell,INeib)
 
 Call MeshPreprocForImplicit(Dim,NC,NFace_Cell,IFace_Cell,INeib,NnonzeroCell,InonzeroCell,InonzeroFace)
 
!Part 9:
 DO cell=1,NC
    CALL PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 END DO
 print*,'find points of cell finished'
 
!Part 10:
 Call GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)

!Part 11:
 !!!Call WallDistance3D(Dim,NC,NFW1,NFW2,FaceType,IDS,X,Y,Z,Xc,Yc,Zc,INW,DW)
 print*,'pre proc finished'
 
!Part 12:
 Call InitMeanFlow3D(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,W0,T0,Mu0,B0,WNP1)

!Part 13:
 Do J=1,NC 
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    W = WNP1(4,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))
    
   !Part 14:
    Temp = GM*P(J)/WNP1(1,J) 

   !Part 15:
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)     

 End Do

!Part 16:
 Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)

!Part 17:
 Ncyc = 0
 Rm   = 10.0
 Naverage = 0
 Time = 0.0
 
!Part 18:
 !!!Call Kw_Init(Dim,NC,MR,WTNP1,Mut)
 !!!Call Ke_Init(Dim,NC,MR,WTNP1,Mut) 
 !!!Call KwSST_Trans_Init3D(Dim,NC,IDS,Minf,Rinf,Wtnp1,Mut,TUinf,Mut0,Gamainf,Kinf,Omegainf,WNP1) 
 !!!Call SA_Init3D(Dim,NC,Mu0,R0,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mut0,Nuhat0,WTNP1,Mut)
 TurbQ(:) =  0.0
 Mut(:)=0.0
 
!Part 19:
 Do While(Rm > ERmx)

   !Part 20:
    Ncyc=Ncyc+1
     
   !Part 21:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)
       WC(5,J) = WNP1(5,J)   
    End Do

   !Part 22:
	Call TimSTP_Turb3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,VOL,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)

   !Part 23:
    !!!DT_min = 1000.0
    !!!Do J=1,NC
    !!!   IF( DT_min > DT(J) ) DT_min = DT(J)
    !!!End Do
    !!!DT(:) = DT_min
    !!!Time = Time + DT_min

   !Part 24:
    Do NS=1,NRKS
   
      !Part 25:
	   RKco=RKJ(NS)

      !Part 26:
       !!!Call FirstOrd_Gradient3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,WNP1,WB,P,GWNP1)
       !!!Call Limiter3D(Dim,NC,NF1,NF2,IDS,GM,XC,YC,ZC,FaceType,WNP1,P,WB,GWNP1,Limit,X,Y,Z,NF)
       !!!Call ConMeanFlow_AUSM_Plus_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,FaceType,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Zc,Limit,Con,X,Y,Z)
       
      !Part 27:
	   !!!Call ConMeanFlow_AUSM_PlusUP3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Minf,WNP1,WB,P,Con)
       Call ConMeanFlow_ScalarDiss_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,Con)
       
      !Part 28:
       Call VelTemp_GradFace3D(Dim,NC,NF,NP,NF1,NF2,NFW1,NFW2,NFS1,NFS2,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,&
                              WNP1,WB,GM,P,NX,NY,NZ,&
                              DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ)
    
      !Part 29:
       Call DifMeanFlowTurbNoWallFu3D(Dim,NC,NF1,NF2,NFW1,NFW2,NF,IDS,GM,PrL,PrT,MR,Mu,Mut,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,&
                                      DVX,DVY,DVZ,DWX,DWY,DWZ,WNP1,WB,TurbQ,Dif)

       !Part 20:
        Res(1,:) = Con(1,:)
        Res(2,:) = Con(2,:) + Dif(2,:)
        Res(3,:) = Con(3,:) + Dif(3,:)
        Res(4,:) = Con(4,:) + Dif(4,:)
        Res(5,:) = Con(5,:) + Dif(5,:)

       !Part 20:
        SOR=0.2

        !!!Call LUSGS(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace,Vol,DT,WNP1&
        !!!         &,WB,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,MR,DQ)

        Call GMRES(Dim,IDS,NC,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace&
                            &,Vol,P,Mu,Mut,PrL,PrT,xc,yc,zc,DT,MR,WNP1,WB,DQ)
  
        !!!Call BLUSGS(Dim,IDS,NC,NF1,NF2,NF,NX,NY,NZ,DA,GM,Res,NnonzeroCell,InonzeroCell,InonzeroFace&
        !!!                    &,Vol,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,DT,MR,WNP1,WB,DQ)
             
        Co = RKco
        WNP1(:,:)=WC(:,:) +Co*sor*DQ(:,:)
        Do J=1,NC

         !Part 21:
	      U = WNP1(2,J)/WNP1(1,J)
          V = WNP1(3,J)/WNP1(1,J)
          W = WNP1(4,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))

         !Part 7:
          Temp = GM*P(J)/WNP1(1,J)

         !Part 8:
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)

        End Do
        
    !!!  !Part 30:
    !!!   Do J=1,NC
    !!!
		  !!!Co = RKco*DT(J)/Vol(J)
    !!!     
    !!!      WNP1(1,J) = WC(1,J) - Co* ( Con(1,J)            )
		  !!!WNP1(2,J) = WC(2,J) - Co* ( Con(2,J) + Dif(2,J) )
    !!!      WNP1(3,J) = WC(3,J) - Co* ( Con(3,J) + Dif(3,J) )
    !!!      WNP1(4,J) = WC(4,J) - Co* ( Con(4,J) + Dif(4,J) )
    !!!      WNP1(5,J) = WC(5,J) - Co* ( Con(5,J) + Dif(5,J) )
    !!!      
    !!!     !Part 31:
	   !!!   U = WNP1(2,J)/WNP1(1,J)
    !!!      V = WNP1(3,J)/WNP1(1,J)
    !!!      W = WNP1(4,J)/WNP1(1,J)
    !!!      P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))
    !!!      
    !!!     !Part 32:
    !!!      Temp = GM*P(J)/WNP1(1,J) 
    !!!      
    !!!     !Part 33:
    !!!      Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)     
    !!!
    !!!   End Do

      !Part 34: 
       Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
       Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)

    End Do !Ns	

   !Part 35:
   !=====================================================================================================================
    !!!Call KWSST_Main3D(Dim,Ncyc,INW,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
    !!!                  NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call KeLB_Main(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
    !!!               NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Wb,WNP1,Mu,P,GM,DUY,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call KwWilcox_Main3D(Dim,NC,NP,INW,IDS,MR,NRKS,RKJ,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,DA,FaceType,&
    !!!                     NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,Z,NX,NY,NZ,XC,YC,ZC,DW,Vol,Mu,DT,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call LES_WALE_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,DW,Mut)
    !!!TurbQ(:) =  0.0

    !!!Call LES_DSmag_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,Mut,TurbQ)
    !!!
    !!!Call KWSST_Transition_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,&
    !!!                             NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,&
    !!!                             DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut,Kinf,Omegainf,Gamainf)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call SA_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
    !!!               NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,WB,WNP1,&
    !!!               Mu0,Mu,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0,WTNP1,Mut)
    !!!TurbQ(:) =  0.0
    
    
    !!!Call KeLB_MainImplicit_f(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,&
    !!!                          NP,IDS,FaceType,INW,XC,YC,ZC,DW,DA,DT,Vol,MR,Wb,WNP1,Mu,P,&
    !!!                          GM,DUY,InonzeroCell,InonzeroFace,NnonzeroCell,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
   !=====================================================================================================================

    !Part 36:
 	 Call ResMass3D(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,'Ncyc:',Ncyc,'Mass eq. Residual:',Rm

    !Part 37:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm

      Call Write_CP3D(Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P)
      Call Write_CF3D(Dim,Minf,Rinf,NFW1,NFW2,X,Y,Z,IDS,DUY,Mu)
      
      Call Write_FlatPlateAnalyticData3D(Dim,NFW1,NFW2,NC,IDS,X,Xc,Yc,Minf,Rinf,WNP1,Mu,DUY)
      
      Call Write_VelocityContour3D(Dim,NC,NP,X,Y,Z,Corn,WNP1)
      
      Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,P,"Pressr.plt")
      Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,Mut,"Muturb.plt")
      
      Call Write_ConservativeVariables3D(Dim,NC,WNP1)
     End If
     
 End Do !Do While
 
!*******************************************************************************************
 End 
!###########################################################################################


