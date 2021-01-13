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
!// This Code designed for analyzing steady Turbulence flow by solving three dimensional Navier-Stocks //!
!// equations. One of the main advantages of this code is that it enjoy several subroutines that       //!
!// assists user to follow code step by step. Mesh and Settings have external files that read through  //!
!// their subroutines in the first step of code. In addition Gauss-Green theorem used for calculating  //!
!// gradients. This code uses wide range wide range of turbulence models ranging for calculating       //!
!// turbulence viscosity Main features of this code is as following:                                   //!
!// 1-Dimension:	                            3D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Turbulent                                              //!
!// 7-Convection Discretization Scheme:	        AUSM+                                                  //!
!//                                             Ausm+Up                                                //!
!//                                             ScalarDiss                                             //!
!// 8-Convection Term Discretization Accuracy:	First order                                            //!
!//                                             Second Order                                           //!
!// 9-Transient Term Discretization Scheme:     Dual Time Stepping                                     //!
!// 10-Steady/Unsteady: 	                    Steady and Unsteady                                    //!
!// 11-Gradient Calculation Scheme:	            Green-Gauss                                            //!
!//                                             Least Square                                           //!
!// 12-Turbulence Model:                        KWSST                                                  //!
!//                                             KeLB                                                   //!
!//                                             LES_WALE                                               //!
!//                                             LES_DSmag                                              //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program 03-DANS_DualTimStp3D
 Implicit None
!===============================
 Integer,Parameter::Dim=500000
!===============================
 Integer::I,J
 Integer::CELL    !  COUNTER for numBER OF CELLS
 Integer::NC !Number of Cells of mesh
 Integer::NP  !Number of Existing Points
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NR  !!Number Of Regions of mesh
 Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NFW1,NFW2 !Index of 1st and last Faces on Wall Boundary 
 Integer::NFF1,NFF2 !Index of 1st and Last Faces on Far-Field Boundary
 Integer::NFI1,NFI2 !Index of 1st and Last Faces on Inflow Boundary 
 Integer::NFS1,NFS2 !Index of 1st ans Last Faces on Symmetry Boundary
 Integer::NFO1,NFO2 !Index of 1st and Last Faces on Inflow Boundary
 Integer::NFIF1,NFIF2 !Index of 1st and last Faces on InterFsce Boundary
 Integer::NRKS !Number of Runge Kutta Stages
 Integer::NS !Number of Steps in rung-kutta method
 Integer::NWrite !Number of Cycle to Write Results
 Integer::Init !Initialize from Available Data in the file(1) or Infinite Flows(0),cell
 Integer::NBP !Number of Boundary Points
 Integer::Unsteady_Moving !!  With this variable first of simulation steady state solution is computed then unsteady moving mesh is computed
 Integer::Test_Case !!Number of test case
 Integer::N_Positions !Number of Positions to Write Results
 Integer::Position_Num !!  Counter of Positions for writing results
 Integer::IStp !index of moving step
 Integer::Stdy_cyc !Number Of Steady Cycles
 Integer::Ncyc !Nymber of Cycles of iterations 
 Integer::II   !!If this term be equal to one steady cycle is completed and unsteady tend to  be started
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::Co    !!  Coefficients of residuals in Rung-kutta steps
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0,V0,W0 !Component of Infinite Flow Velocity
 Real(8)::E0 !Internal Energy of Infinite Flow
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::TT !Total Temprature
 Real(8)::B0 !Sauterland Constant
 Real(8)::RM !Residual of Mass Equation
 Real(8)::MR !Much Number over Reynolds Number of infinite Flow
 Real(8)::PrL !Prantle Number for Laminar Flows
 Real(8)::PrT !Prantle Number for Turbulent Flows
 Real(8)::T0 !Temperature of Infinite Flow
 Real(8)::MU0  !Molecular Viscosity of Infinite Flow
 Real(8)::Mut0 !Eddy Viscosity of Infinite Flow
 Real(8)::U    !  Velocity in x Direction
 Real(8)::V    !  Velocity in Y Direction
 Real(8)::RKco !  Rung-kutta Coefficients
 Real(8)::W   !  Velocity in Z Direction
 Real(8)::Temp  !  Temprature  
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::DT_REAL !Minimum Time Step*Time_Coe
 Real(8)::Real_Time !  Real Time Passed For converging
 Real(8)::Total_Time !  Total time Considerd For Solution By User
 Real(8)::CL !Lift Coefficient
 Real(8)::CD !Drag Coefficient
 Real(8)::phiWrite  !!  Local Varible for Counter of Positions for writing results
 Real(8)::Delphi !!  Differentiation In Phi Direction
 Real(8)::phi    !  Momentary Angle Of Phase
 Real(8)::phi1  !  Represent Residual of Division phi per 360
 Real(8)::Omega  !!Frequency of Flapping Airfoil
 Real(8)::PI  !!Mathematical constant that equals to 3.14159
 Real(8)::X1 !!First Point Coordinate of First Line
 Real(8)::H  !  Local Variable that Constitute Horizental Axis of Lift Coefficient Graph
 Real(8)::StdyERmx !!Value Of ERror considered for mass equation by user
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8)::UnStdyResMas !!Residual of Real Time Iteration
 Real(8)::StdyResMas !!  Residual of Imaginary Iteration
 Real(8)::AA ! Local Variable for Area of each cell
 Real(8)::PlanForm_Area
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Integer,Dimension(1:Dim)::NFace_Cell !Number of Faces Forming a Cell
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:6,1:Dim)::IFace_Cell !Index of Faces Forming a Cell
 Integer,Dimension(1:8,1:Dim)::Corn !Corners point index of Each Cell 
 Real(8),Dimension(1:40)::phi_Write !!  Phase of Positions to Write Results
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xc,Yc,Zc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::Vol !Volume of each Cell
 Real(8),Dimension(1:Dim)::NX,NY,NZ !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::DT !Time step
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:3,1:Dim)::FACE_VELOCITY !Velocity of Faces
 Real(8),Dimension(1:5,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WN !Conservative Values at (N)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WNM1 !Conservative Values at (N-1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(1:5,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:5,1:Dim)::Dif ! Diffusion Term of Mean flow Equations
 Real(8),Dimension(1:5,1:Dim)::Res !Residual of Equations
 Real(8),Dimension(1:5,1:Dim)::WTNP1 !Turbulence Variables at new time step
 Real(8),Dimension(1:5,1:Dim)::WTN !Turbulent Conservative Values at (N)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WTNM1 !Turbulent Conservative Values at (N-1)th Time Step
 Real(8),Dimension(1:6,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Integer::Naverage   !Number of Averaging Steps=0
 Integer,Dimension(1:Dim)::INW  !Index of Nearest Wall
 Real(8),Dimension(1:Dim)::DW  !Distance to Nearest Wall
 Real(8),Dimension(1:Dim)::MU !Molecular Viscosity
 Real(8),Dimension(1:Dim)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ !Derivative component of Velocity and Temperature
 Real(8),Dimension(1:Dim)::TurbQ !Turbulence Quantity (i.e. k in turbulence model)
 Real(8),Dimension(1:Dim)::Mut  !Turbulence Viscosity
 Real(8),allocatable,Dimension(:)::sum_umean !Summation of Averegaed X Component Velocity (u)
 Real(8),allocatable,Dimension(:)::sum_vmean  !Summation of Averegaed Y Component Velocity (v)
 Real(8),allocatable,Dimension(:)::sum_uvmean !Summation of Averegaed u*v Velocities
 Real(8),allocatable,Dimension(:)::sum_u2mean !Summation of Averegaed u*u Velocities
 Real(8),allocatable,Dimension(:)::sum_v2mean !Summation of Averegaed V*V Velocities
!***************************************** Main ********************************************
 PI=4.*Atan(1.)
!Part 1:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3: 

!Part 4:
 Call MeshBC3D(Dim,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 6:
 Call FaceOfCell(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 Do Cell=1,NC
    Call PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 End Do
 
!Part 7:
 Call GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)

 
 Call WallDistance3D(Dim,NC,NFW1,NFW2,FaceType,IDS,X,Y,Z,Xc,Yc,Zc,INW,DW)
 
!Part 8:
 Call InitMeanFlow3D(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,W0,T0,Mu0,B0,WNP1)
 
!Part 9:
 Do J=1,NC 
    U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    W = WNP1(4,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))

   !Part 7:
    Temp = GM*P(J)/WNP1(1,J)

   !Part 8:
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)
    
    WN(1,J) = WNP1(1,J)
    WN(2,J) = WNP1(2,J)
    WN(3,J) = WNP1(3,J)
    WN(4,J) = WNP1(4,J) 
    WN(5,J) = WNP1(5,J)
     
     WTN(1,J) = WTNP1(1,J)
     WTN(2,J) = WTNP1(2,J)
 End Do
 
!Part 11:
 Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 
 !Part 9:
 Call Ke_Init(Dim,NC,MR,WTNP1,Mut)  
 !!!Call Kw_Init(Dim,NC,MR,WTNP1,Mut) 
 
 !Part 10:
 Ncyc = 0
 Real_Time = 0.0

!Part 13:
 Do While(Real_Time < TotTime)
     
   !Part 14:
    Ncyc=Ncyc+1
     
   !Part 15:
    Do J=1,NC
        WNM1(1,J) = WN(1,J)
        WNM1(2,J) = WN(2,J)
        WNM1(3,J) = WN(3,J)
        WNM1(4,J) = WN(4,J)
        WNM1(5,J) = WN(5,J)
        
        WN(1,J) = WNP1(1,J)
        WN(2,J) = WNP1(2,J)
        WN(3,J) = WNP1(3,J)
        WN(4,J) = WNP1(4,J)
        WN(5,J) = WNP1(5,J)
         
        WTNM1(1,J) = WTN(1,J)
        WTNM1(2,J) = WTN(2,J)
         
        WTN(1,J) = WTNP1(1,J)
        WTN(2,J) = WTNP1(2,J)
    End Do
    
   !Part 16:
    Call TimSTP_Turb3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,VOL,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)
    
   !Part 17:
    DT_Real = 1000.0
    Do J=1,NC
        If( DT_Real > DT(J) ) DT_Real = DT(J)
    End Do
    DT_Real = DT_Real * Time_Coe
    
   !Part 18:
    Real_Time = Real_Time+DT_Real
    
   !Part 26:
    RM = 100.0
    Stdy_cyc = 0
    
   !Part 27:
    Do While(RM > ERmx .or. Stdy_cyc<2000)
        
       !Part 28:
        Stdy_cyc = Stdy_cyc + 1
        
       !Part 29:
        Do J=1,NC
            WC(1,J) = WNP1(1,J)
            WC(2,J) = WNP1(2,J)
            WC(3,J) = WNP1(3,J)
            WC(4,J) = WNP1(4,J)
            WC(5,J) = WNP1(5,J)
        End Do
        
       !Part 30:
        Do NS=1,NRKS
            
           !Part 31:
            RKco=RKJ(NS)
            
           !Part 32:
            Call ConMeanFlow_ScalarDiss_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,Con)

            Call VelTemp_GradFace3D(Dim,NC,NF,NP,NF1,NF2,NFW1,NFW2,NFS1,NFS2,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,&
                                    WNP1,WB,GM,P,NX,NY,NZ,&
                                    DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ)

            Call DifMeanFlowTurbNoWallFu3D(Dim,NC,NF1,NF2,NFW1,NFW2,NF,IDS,GM,PrL,PrT,MR,Mu,Mut,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,&
                                           DVX,DVY,DVZ,DWX,DWY,DWZ,WNP1,WB,TurbQ,Dif)

           !Part 33:
            Do J=1,NC
                
                AA = 1./Vol(J)
                Res(1,J) = -( Con(1,J)            )*AA
                Res(2,J) = -( Con(2,J) + Dif(2,J) )*AA
                Res(3,J) = -( Con(3,J) + Dif(3,J) )*AA
                Res(4,J) = -( Con(4,J) + Dif(4,J) )*AA
                Res(5,J) = -( Con(5,J) + Dif(5,J) )*AA
                
                Temp = 2*DT_Real+3*RKco*DT(J)
                WNP1(1,J) = ( 2*WC(1,J)*DT_Real + RKco*DT(J)*(4*WN(1,J) - WNM1(1,J) + 2*Res(1,J)*DT_Real) ) / Temp
                WNP1(2,J) = ( 2*WC(2,J)*DT_Real + RKco*DT(J)*(4*WN(2,J) - WNM1(2,J) + 2*Res(2,J)*DT_Real) ) / Temp
                WNP1(3,J) = ( 2*WC(3,J)*DT_Real + RKco*DT(J)*(4*WN(3,J) - WNM1(3,J) + 2*Res(3,J)*DT_Real) ) / Temp
                WNP1(4,J) = ( 2*WC(4,J)*DT_Real + RKco*DT(J)*(4*WN(4,J) - WNM1(4,J) + 2*Res(4,J)*DT_Real) ) / Temp	
                WNP1(5,J) = ( 2*WC(5,J)*DT_Real + RKco*DT(J)*(4*WN(5,J) - WNM1(5,J) + 2*Res(5,J)*DT_Real) ) / Temp
                
                U = WNP1(2,J)/WNP1(1,J)
                V = WNP1(3,J)/WNP1(1,J)
                W = WNP1(4,J)/WNP1(1,J)
                P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))

               !Part 7:
                Temp = GM*P(J)/WNP1(1,J)

               !Part 8:
                Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)
       
            End Do !J
            
           !Part 34: 
            Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
            Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
            Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
            Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
            Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
            
        End Do !Ns
        
       !Part 35:
 	    Call ResMass3D(Dim,NC,WNP1,WC,DT,Stdy_cyc,RM)
	    print*,'Stdy_cyc: ',Stdy_cyc,RM
      
    End Do !Steady

    !!!Call KeLB_Main3D_DualTimStp(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
    !!!                            NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Wb,WNP1,Mu,P,GM,DUY,WTNM1,WTN,DT_Real,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    Call KWSST_Main3D_DualTimStp(Dim,Ncyc,INW,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                      NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,WTNM1,WTN,DT_Real,Mu,WTNP1,Mut)
    TurbQ(:) =  WTNP1(1,:)
    
    !!!Call LES_WALE_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,DW,Mut)
    !!!TurbQ(:) =  0.0

    !!!Call LES_DSmag_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,Mut,TurbQ)
    
    
    Call CalculateLESVariables3D(Dim,Init,Ncyc,NC,WNP1,Real_Time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean)

   !Part 37:
    Print*,Ncyc,DT_Real,Real_Time
         
   !Part 36:
    If( Mod(Ncyc,NWrite)==0 )Then
         
    !Part 38:
    Call Write_CF3DUnsteady(Dim,NFW1,NFW2,IDS,FaceType,X,Y,Z,Minf,Rinf,Mu,DUY,Naverage,SumCF)
    Call Write_CP3DUnsteady(Dim,NFW1,NFW2,IDS,FaceType,X,Y,Z,Minf,P,GM,Naverage,SumCP)

    Call Write_CP3D(Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P)
    Call Write_CF3D(Dim,Minf,Rinf,NFW1,NFW2,X,Y,Z,IDS,DUY,Mu)

    Call Write_VelocityContour3D(Dim,NC,NP,X,Y,Z,Corn,WNP1)

    Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,P,"Pressr.plt")
    Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,Mut,"Muturb.plt")

    Call Write_LESVariables3D(Dim,NC,NP,Corn,X,Y,Z,Real_Time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean)

    Call Write_ConservativeVariables3D(Dim,NC,WNP1)
         
    !Part 39:
    Call PressLiftDragCo_3D(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,NZ,DA,Unsteady_Moving,PlanForm_Area,CL,CD)
    Write(10,*)  X1,CL
    Write(100,*) X1,CD
         
    End If
     

 End Do !Do While
!*********************************************************************************************
 End 
!###########################################################################################


