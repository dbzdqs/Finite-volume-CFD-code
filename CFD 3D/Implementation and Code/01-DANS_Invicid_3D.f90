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
!// This Code designed for solving Inviscid form of two dimensional simplified Navier-Stocks equations //!
!// known as Euler equations. One of the main advantages of this code is that it enjoy several         //!
!// subroutinesthat assists user to follow code step by step. Mesh and Settings have external files    //!
!// that read through their subroutines in the first step of code. User can test and compare various   //!
!// methods for convective terms discretization. For example code can switch from central methods of   //!
!// discretizing to upwind methods only by selecting and replacing proper subroutine. Main features of //!
!// this code is as following:                                                                         //!
!// 1-Dimension:	                            3D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Invicid                                                //!
!// 7-Convection Discretization Scheme:	        AUSM+                                                  //!
!//                                             Ausm+Up                                                //!
!//                                             ScalarDiss                                             //!
!// 8-Convection Term Discretization Accuracy:	First order                                            //!
!//                                             Second order                                           //!
!// 9-Transient Term Discretization Scheme:     Explicit(Rung-Kutta)                                   //!
!// 10-Steady/Unsteady: 	                    Steady and Unsteady                                    //!
!// 11-Gradient Calculation Scheme:	            Green-Gauss                                            //!
!//                                             Least Squar                                            //!
!// 12-Turbulence Model:	                    non                                                    //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Invicid_3D
 Implicit None
!===============================
 Integer,Parameter::Dim=200000
!===============================
 
 Integer::I,J
 Integer::NS    !  Number of Steps in rung-kutta method
 Real(8)::U    !  Velocity in x Direction
 Real(8)::V    !  Velocity in Y Direction
 Real(8)::W    !  Velocity in Z Direction
 Real(8)::RKco !  Rung-kutta Coefficients
 Real(8)::co   !  Coefficients of residuals in Rung-kutta steps
 Integer::Cell    !  COUNTER for number of Cells
 Integer::NC !Number of Cells of mesh
 Integer::NP  !Number of Existing Points
 Integer::NF !Number of Faces Constructing Mesh
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
 Real(8)::RM !Residual of Mass Equation 
 Integer::Naverage   !Number of Averaging Steps
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0,V0,W0   !Component of Infinite Flow Velocity
 Real(8)::E0 !Internal Energy of Infinite Flow 
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TT !Total Temprature
 Real(8)::time       !Simulation Time
 Real(8)::DT_min !Minimum Time Step !Minimum Time Step
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,allocatable,Dimension(:,:)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),allocatable,Dimension(:,:)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),allocatable,Dimension(:,:)::WC !Constant values Of Rung-kutta Method
 Real(8),allocatable,Dimension(:,:)::Con !Convection Term of Mean flow Equations
 Real(8),allocatable,Dimension(:)::X,Y,Z !Coordinate of Points Constructing Mesh 
 Real(8),allocatable,Dimension(:)::Xc,Yc,Zc !Coordinate of Center of Element
 Real(8),allocatable,Dimension(:)::Vol !Volume of each Cell !Volume of each Cell
 Real(8),allocatable,Dimension(:)::NX,NY,NZ !Normal Vectors of each Face
 Real(8),allocatable,Dimension(:)::DA   !Area of each Face
 Real(8),allocatable,Dimension(:)::DT !Time step
 Real(8),allocatable,Dimension(:)::P !Pressure
 Real(8),allocatable,Dimension(:,:)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,allocatable,Dimension(:)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,allocatable,Dimension(:)::N_Corn !Number of point Forming a Cell
 Integer,allocatable,Dimension(:,:)::IFace_Cell !Index of Faces Forming a Cell
 Integer,allocatable,Dimension(:)::NFace_Cell !Number of Faces Forming a Cell 
 Integer,allocatable,Dimension(:,:)::Corn !Corners point index of Each Cell 
 Real(8),allocatable,Dimension(:)::sum_umean !Summation of Averegaed X Component Velocity (u)
 Real(8),allocatable,Dimension(:)::sum_vmean  !Summation of Averegaed Y Component Velocity (v)
 Real(8),allocatable,Dimension(:)::sum_uvmean !Summation of Averegaed u*v Velocities
 Real(8),allocatable,Dimension(:)::sum_u2mean !Summation of Averegaed u*u Velocities
 Real(8),allocatable,Dimension(:)::sum_v2mean !Summation of Averegaed V*V Velocities
 Real(8),allocatable,Dimension(:,:)::Limit !Limiter for high order calculation 
 Real(8),allocatable,Dimension(:,:,:)::GWNP1 !Gradients of Conservative Values at (N+1)th Time Step
!***************************************** Main ********************************************
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
 allocate( WB(1:6,1:Dim)         )
 allocate( FaceType(1:Dim)       )
 allocate( N_Corn(1:Dim)         )
 allocate( IFace_Cell(1:6,1:Dim) )
 allocate( NFace_Cell(1:Dim)     )
 allocate( Corn(1:8,1:Dim)       )
 allocate( DT(1:Dim)             )
 allocate( P(1:Dim)              )

 !!!allocate( Limit(1:5,1:Dim)      )
 !!!allocate( GWNP1(1:3,1:5,1:Dim)  )

!Part 1:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)

!Part 3:
 Call MeshBC3D(Dim,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 4:
 Call FaceOfCell(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 
!Part 6:
 DO cell=1,NC
    CALL PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 END DO

!Part 7:
 Call GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)

!Part 8:
 Call InitMeanFlow_Inviscid3D(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,W0,WNP1)
  
!Part 9: 
 Do J=1,NC
   
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    W = WNP1(4,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))

 End Do

!Part 11:
 Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)

!Part 12:
 Ncyc = 0
 Rm   = 10.0

!Part 13:
 Do While(Rm > ERmx)

   !Part 14:
    Ncyc=Ncyc+1
     
   !Part 15:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)
       WC(5,J) = WNP1(5,J)   
    End Do

   !Part 16:
	Call TimSTP_Inviscid3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,DA,CFLx,GM,P,WNP1,WB,DT)

   !Part 17:
    Do NS=1,NRKS
   
      !Part 18:
	   RKco=RKJ(NS)

       !Part 19:
       !!!Call Limiter3D(Dim,NC,NF1,NF2,IDS,GM,XC,YC,ZC,FaceType,WNP1,P,WB,GWNP1,Limit,X,Y,Z,NF)
       !!!Call ConMeanFlow_AUSM_Plus_HO(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,FaceType,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Zc,Limit,Con,X,Y,Z)
	   Call ConMeanFlow_AUSM_PlusUP3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Minf,WNP1,WB,P,Con)
       !!!Call ConMeanFlow_ScalarDiss_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,Con)
       
      !Part 20:
       Do J=1,NC

		  Co = RKco*DT(J)/Vol(J)
         
          WNP1(1,J) = WC(1,J) - Co* Con(1,J)
		  WNP1(2,J) = WC(2,J) - Co* Con(2,J)
          WNP1(3,J) = WC(3,J) - Co* Con(3,J)
          WNP1(4,J) = WC(4,J) - Co* Con(4,J)
          WNP1(5,J) = WC(5,J) - Co* Con(5,J)

         !Part 21:
	      U = WNP1(2,J)/WNP1(1,J)
          V = WNP1(3,J)/WNP1(1,J)
          W = WNP1(4,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))
    
       End Do

      !Part 23: 
       Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
       Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
       
    End Do !Ns	
    
    !Part 24:
 	 Call ResMass3D(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,Ncyc,Rm

    !Part 25:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm

      Call Write_CP3D(Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P)
      
      Call Write_VelocityContour3D(Dim,NC,NP,X,Y,Z,Corn,WNP1)
      
      Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,P,"Pressr.plt")
      
      Call Write_ConservativeVariables3D(Dim,NC,WNP1)
     End If
     
 End Do !Do While
 
!*********************************************************************************************
 End 
!###########################################################################################


