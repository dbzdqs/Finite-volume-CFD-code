!DDDDDDDDDDDDDDDDDDDDDDDDDDDAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNSSSSSSSSSSSSSSSSSSSSSSSSSS
!//             /////////////       ////////////////    ////////     //////    ////////////////        //!
!//             /////////////      ////////////////    //////////   //////    ////////////////         //!
!//            /////    /////     //////    //////    //////////// //////    /////                     //!
!//           /////    //////    ////////////////    ///////////////////    ////////////////           //!
!//          /////    //////    ////////////////    ////// ////////////               /////            //!
!//         ///////////////    //////    //////    //////   //////////    ////////////////             //!
!//       ///////////////     //////    //////    //////     ////////    ////////////////              //!
!//          Developer            Assistant    in      Numerical             Sciences                  //!
!//----------------------------------------------------------------------------------------------------//!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace Department, Amirkabir University of Technology           //!
!// Chief Developer: N. msnkre, Aerospace eng., Amirkabir University of Technology                     //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.DANS.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!// This code are designed to simulate the effect of DBD plasma actuator of flow field. The body force //! 
!// of this kinds of actuator is a body force which is coupled with N.S. equations as a source term.   //!
!// The DBD actuator is simulated by Shyy model. The flow solver for solving N.S. equation is a finite //!
!// volume flow solver code using unstructured mesh. There many turbulence model and discretization    //!
!// scheme for numerical modeling of flow field.                                                       //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Turb
 Implicit None
!===============================
 Integer,Parameter::Dim=120000
!===============================
 Integer::I,J,NS  
 Real(8)::Co,U,V,RKco,Temp
          
 Integer::NC !Number of Cells of mesh
 Integer::NP  !Number of Existing Points
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NFW1,NFW2 !Index of 1st and last Faces on Wall Boundary 
 Integer::NFF1,NFF2 !Index of 1st and Last Faces on Far-Field Boundary
 Integer::NFI1,NFI2 !Index of 1st and Last Faces on Inflow Boundary 
 Integer::NFS1,NFS2 !Index of 1st ans Last Faces on Symmetry Boundary
 Integer::NFO1,NFO2 !Index of 1st and Last Faces on Inflow Boundary
 Integer::NFIF1,NFIF2 !Index of 1st and last Faces on InterFsce Boundary
 Integer::NRKS !Number of Runge Kutta Stages
 Integer::NWrite   !Number of Cycle to Write Results
 Integer::NR   !Number Of Regions of mesh
 Integer::Init !Initialize from Available Data in the file(1) or Infinite Flows(0)
 Integer::Ncyc !Nymber of Cycles of iterations 
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0   !Sound Speed of Infinite Flow
 Real(8)::U0   !Infinite Flow Velocity in X Direction
 Real(8)::V0   !Infinite Flow Velocity in Y Direction
 Real(8)::E0 !Internal Energy of Infinite Flow
 Real(8)::T0
 Real(8)::H0 !Enthalpy of Infinite Flow
 Real(8)::B0 !Sauterland Constant
 Real(8)::MU0  !Molecular Viscosity of Infinite Flow
 Real(8)::Mut0 !Eddy Viscosity of Infinite Flow 
 Real(8)::Rinf !Reynolds Number of infinite Flow  
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::TT !Total Temprature 
 Real(8)::MR !Much Number over Reynolds Number of infinite Flow 
 Real(8)::PrL !Prantle Number for Laminar Flows 
 Real(8)::PrT !Prantle Number for Turbulent Flows 
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::RM !Residual of Mass Equation
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(1:4,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:4,1:Dim)::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(1:4,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:4,1:Dim)::Dif ! Diffusion Term of Mean flow Equations
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xc,Yc !Coordinate of Center of Element 
 Real(8),Dimension(1:Dim)::A !Area of each cell
 Real(8),Dimension(1:Dim)::NX,NY !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::DW  !Distance to Nearest Wall
 Real(8),Dimension(1:Dim)::DT !Time step
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:Dim)::MU !Molecular Viscosity
 Real(8),Dimension(1:Dim)::Mut  !Turbulence Viscosity
 Real(8),Dimension(1:Dim)::DUX,DUY,DVX,DVY,DTX,DTY !Derivative component of Velocity and Temperature
 Real(8),Dimension(1:5,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,Dimension(1:Dim)::INW  !Index of Nearest Wall
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Integer,Dimension(1:4,1:Dim)::Corn !Corners point index of Each Cell 
 Real(8),Dimension(1:Dim)::TurbQ !Turbulence Quantity (i.e. k in turbulence model)
 Real(8),Dimension(1:2,1:Dim)::WTNP1 !Turbulence Variables at new time step

 Real(8)::freq
 Real(8)::PlasmaHeight      !height of Plasma in Y-Axis (ND)
 Real(8)::Plasmawidth       !width of Plasma in X-Axis (ND)
 Real(8)::appliedVoltage    !applied voltage (kv)
 Real(8)::PlasmaGap         !Distance Between the plates (cm)
 Real(8)::Roh_c		        !Density of Electron (/cm^2)
 Real(8)::e_c		        !charge of Electron (c)
 Real(8)::Eb			    !Breakdown Electric Field Strength (kv/cm)
 Real(8)::DischargeTime     !Discharge Time (sec)
 
 Real(8)::roh_c_max
 Real(8)::phi_max
 Real(8)::f_t
 Real(8)::f_prim_t
 Real(8)::sigma_L
 
 Integer,Dimension(1:Dim)::delta
 Real(8),Dimension(1:Dim)::F_DBD_x,F_DBD_y
 
 Integer::Naverage   !Number of Averaging Steps=0
 Real(8),Dimension(1:Dim)::sum_umean=0.0,sum_vmean=0.0,umean=0.0,vmean=0.0,SumCP=0.0,SumCF=0.0
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
 Call PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,IDS,Corn)
 
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
 !!!Call Ke_Init(Dim,NC,MR,WTNP1,Mut)  
 Call Kw_Init(Dim,NC,MR,WTNP1,Mut) 
 
!Part 10:
 Call PlasmaAfectedCell_General(Dim,NC,Xc,Yc,Delta)
 Call PlasmaShayyParameters(freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,Roh_c,e_c,Eb,DischargeTime)

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
	   Call DifMeanFlow_TurbNoWallFu(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,GM,PrL,PrT,NX,NY,MR,Mu,Mut,&
                                WNP1,WTNP1(1,:),WB,DUX,DUY,DVX,DVY,DTX,DTY,Dif)

      !Part 21:
       Do J=1,NC

		  Co = RKco*DT(J)/A(J)

          WNP1(1,J) = WC(1,J) - Co*( Con(1,J)          )
		  WNP1(2,J) = WC(2,J) - Co*( Con(2,J)+Dif(2,J) - F_DBD_x(J) )
          WNP1(3,J) = WC(3,J) - Co*( Con(3,J)+Dif(3,J) - F_DBD_y(J) )
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
     
     Call KWSST_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                     NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
     
     !!!Call KWSST_Sust_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
     !!!                NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
    
     !!!Call KWSST_V_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
     !!!                NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,Mu0,Wb,WNP1,Mu,WTNP1,Mut)    
    
     
    !Part 26:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,Ncyc,Rm

     Call UnsteadyVar2D(Dim,WNP1,Naverage,sum_umean,sum_vmean,umean,vmean)
     
    !Part 27:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
      !!!Call Write_CP3DUnsteady(Dim,NFW1,NFW2,IDS,X,Y,Minf,P,GM,Naverage,SumCP)
      Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
      !!!Call Write_CF2DUnsteady(Dim,NFW1,NFW2,IDS,X,Y,Minf,Rinf,Mu,DUY,Naverage,SumCF)
      Call Write_CF(Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu)
      Call Write_ConservativeVariables(Dim,NC,WNP1)
	  Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mu,"MolVis.plt")
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mut,"TubVis.plt")
	 End If

 End Do !Do While
!*********************************************************************************************
 End 
!###########################################################################################


