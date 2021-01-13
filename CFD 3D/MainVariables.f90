Module MainVariables


 Integer::nLocalArrayFace
 Integer::nLocalArrayCell
 Real(8),Dimension(:,:),pointer::RealLocalArrayFace
 Real(8),Dimension(:,:),pointer::RealLocalArrayCell
 Integer,Dimension(:,:),pointer::IntegerLocalArrayFace
 Integer,Dimension(:,:),pointer::IntegerLocalArrayCell
 
 
 Integer::NS      !  Number of Steps in rung-kutta method
 Real(8)::U       !  Velocity in x Direction
 Real(8)::V       !  Velocity in Y Direction
 Real(8)::W       !  Velocity in Z Direction
 Real(8)::RKco    !  Rung-kutta Coefficients
 Real(8)::Temp    !  Temprature
 Real(8)::co      !  Coefficients of residuals in Rung-kutta steps
 Integer::CELL    !  COUNTER for num 
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
 Integer::NR   !Number Of Regions of mesh
 Integer::NRKS !Number of Runge Kutta Stages
 Integer::NWrite   !Number of Cycle to Write Results
 Integer::Init !Initialize from Available Data in the file(1) or Infinite Flows(0)
 Integer::Ncyc !Nymber of Cycles of iterations 
 Integer::Naverage   !Number of Averaging Steps
 Real(8)::GM   !Gama Constant (Specific Heat Ratio),Co
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0,V0,W0   !Component of Infinite Flow Velocity
 Real(8)::E0 !Internal Energy of Infinite Flow
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::TT !Total Temprature
 Real(8)::B0 !Sauterland Constant
 Real(8)::RM !Residual of Mass Equation
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::MR !Much Number over Reynolds Number of infinite Flow
 Real(8)::PrL !Prantle Number for Laminar Flows
 Real(8)::PrT !Prantle Number for Turbulent Flows
 Real(8)::T0 !Temperature of Infinite Flow
 Real(8)::MU0  !Molecular Viscosity of Infinite Flow
 Real(8)::time       !Simulation Time
 Real(8)::Mut0 !Eddy Viscosity of Infinite Flow
 Real(8)::DT_min !Minimum Time Step
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::Kinf ! Infinite Value of Kinetic Energy
 Real(8)::Omegainf ! Infinite Value of specific rate of dissipation
 Real(8)::Gamainf ! Infinite Value of inttermetency
 Real(8)::Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1 !Spalart-Allmaras Parameters
 Real(8)::Nuhat0  !Turbulence Variable in Spalart-Allmaras
 Real(8)::TUinf ! TUrbulence Intensity of Infinite flow
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient

 Integer,Dimension(:),pointer::NFR !Number of  Face of each Regions
 Integer,Dimension(:),pointer::BC  !Boundary Condition index  
 Integer,Dimension(:,:),pointer::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(:,:),pointer::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(:,:),pointer::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(:,:),pointer::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(:,:),pointer::Dif ! Diffusion Term of Mean flow Equations
 Real(8),Dimension(:),pointer::X,Y,Z !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(:),pointer::Xc,Yc,Zc !Coordinate of Center of Element
 Real(8),Dimension(:),pointer::Vol !Volume of each Cell
 Real(8),Dimension(:),pointer::NX,NY,NZ !Normal Vectors of each Face
 Real(8),Dimension(:),pointer::DA   !Area of each Face
 Real(8),Dimension(:),pointer::DT !Time step
 Real(8),Dimension(:),pointer::P !Pressure
 Real(8),Dimension(:),pointer::MU !Molecular Viscosity
 Real(8),Dimension(:),pointer::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ !Derivative component of Velocity and Temperature
 Real(8),Dimension(:),pointer::TurbQ !Turbulence Quantity (i.e. k in turbulence model)
 Real(8),Dimension(:,:),pointer::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,Dimension(:),pointer::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(:),pointer::N_Corn !Number of point Forming a Cell
 Integer,Dimension(:,:),pointer::IFace_Cell !Index of Faces Forming a Cell
 Integer,Dimension(:),pointer::NFace_Cell !Number of Faces Forming a Cell
 Integer,Dimension(:,:),pointer::Corn !Corners point index of Each Cell 
 Integer,Dimension(:),pointer::INW  !Index of Nearest Wall
 Real(8),Dimension(:),pointer::DW  !Distance to Nearest Wall
 Real(8),Dimension(:),pointer::Mut  !Turbulence Viscosity
 Real(8),Dimension(:,:),pointer::WTNP1 !Turbulence Variables at new time step
 Real(8),Dimension(:,:),pointer::WTB !Turbulence Variable at Boundary Faces
 Real(8),Dimension(:),pointer::sumCP !
 Real(8),Dimension(:),pointer::sumCF !
 Real(8),Dimension(:),pointer::sum_umean !Summation of Averegaed X Component Velocity (u)
 Real(8),Dimension(:),pointer::sum_vmean  !Summation of Averegaed Y Component Velocity (v)
 Real(8),Dimension(:),pointer::sum_uvmean !Summation of Averegaed u*v Velocities
 Real(8),Dimension(:),pointer::sum_u2mean !Summation of Averegaed u*u Velocities
 Real(8),Dimension(:),pointer::sum_v2mean !Summation of Averegaed V*V Velocities
 
 Real(8),Dimension(:,:),pointer::Limit !Limiter for high order calculation                        !(1:5,1:Dim)
 Real(8),Dimension(:,:,:),pointer::GWNP1 !Gradients of Conservative Values at (N+1)th Time Step   !(1:3,1:5,1:Dim)

 
 
 
    
    
 contains
    
    
 subroutine MeshVarConstructir(FileName)

 Character *8 FileName
 Open(1,File=FileName)
 Read(1,*) MeshDim
 Read(1,*) NP   
 Read(1,*) NC
 Read(1,*) NF
 Read(1,*) NR
 Close(1)
 
 allocate( IntegerLocalArrayFace(1:nLocalArrayFace,1:NF) ) 
 allocate( RealLocalArrayFace(1:nLocalArrayFace,1:NF) ) 
 allocate( IntegerLocalArrayCell(1:nLocalArrayCell,1:NF) ) 
 allocate( RealLocalArrayCell(1:nLocalArrayCell,1:NF) ) 

 allocate( NFR(1:NR)            )
 allocate( BC(1:NR)             )
 allocate( IDS(1:6,1:NF)        )
 allocate( WNP1(1:5,1:NC)       )
 allocate( WC(1:5,1:NC)         )
 allocate( Con(1:5,1:NC)        )
 allocate( X(1:NP)              )
 allocate( Y(1:NP)              )
 allocate( Z(1:NP)              )
 allocate( XC(1:NC)             )
 allocate( YC(1:NC)             )
 allocate( ZC(1:NC)             )
 allocate( Vol(1:NC)            )
 allocate( NX(1:NF)             )
 allocate( NY(1:NF)             )
 allocate( NZ(1:NF)             )
 allocate( DA(1:NF)             )
 allocate( TurbQ(1:NC)          )
 allocate( WB(1:6,1:NF)         )
 allocate( FaceType(1:NF)       )
 allocate( N_Corn(1:NC)         )
 allocate( IFace_Cell(1:6,1:NC) )
 allocate( NFace_Cell(1:NC)     )
 allocate( Corn(1:8,1:NC)       )
 allocate( DT(1:NC)             )
 allocate( P(1:NC)              )
 
!Part 2:
 allocate( INW(1:NC)       )
 allocate( DW(1:NC)        ) 
 allocate( Dif(1:5,1:NC)   )
 allocate( Mu(1:NC)        )
 allocate( Mut(1:NC)       )
 allocate( DUX(1:NC)       )
 allocate( DUY(1:NC)       )
 allocate( DUZ(1:NC)       )
 allocate( DVX(1:NC)       )
 allocate( DVY(1:NC)       )
 allocate( DVZ(1:NC)       )
 allocate( DWX(1:NC)       )
 allocate( DWY(1:NC)       )
 allocate( DWZ(1:NC)       )
 allocate( DTX(1:NC)       )
 allocate( DTY(1:NC)       )
 allocate( DTZ(1:NC)       )
 allocate( WTNP1(1:2,1:NC) )
 allocate( WTB(1:2,1:NC)   )
 
!Part 3:
 allocate( SumCF(1:NC)      )
 allocate( SumCP(1:NC)      )
 allocate( sum_umean(1:NC)  )
 allocate( sum_vmean(1:NC)  )
 allocate( sum_uvmean(1:NC) )
 allocate( sum_u2mean(1:NC) )
 allocate( sum_v2mean(1:NC) )

 end 
 
 
End Module MainVariables