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
!// This Code designed for solving unsteady inviscid form of two dimensional Euler equations. One of   //!
!// the main advantages of this code is that it enjoy several subroutines that assists user to follow  //!
!// code step by step. Governing equations are formulated based on Arbitary Lagrangian-Eulerian (ALE)  //!
!// approach. Mesh and Settings have external files that read through their subroutines in the first   //!
!// step of code. In addition Dual time stepping scheme used for time marching .Furthermore Radial     //!
!// Base Function (RBF) method used for moving mesh. Main features of this code is as following:       //!
!// 1-Dimension:	                            3D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Invicid/ALE                                            //!
!// 7-Convection Discretization Scheme:	        AUSM                                                   //!
!//                                             AUSM+                                                  //!
!//                                             Ausm+Up                                                //!
!//                                             ScalarDiss                                             //!
!// 8-Convection Term Discretization Accuracy:	First order                                            //!
!// 9-Transient Term Discretization Scheme:     Dual Time Stepping                                     //!
!// 10-Steady/Unsteady: 	                    Unsteady                                               //!
!// 11-Gradient Calculation Scheme:	            non                                                    //!
!// 12-Turbulence Model:	                    non                                                    //!
!// 13-Moving Mesh Method:	                    RBF/Linear Spring                                      //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_ALEInviscid3D_DualTime
 Implicit None
!===============================
 Integer,Parameter::Dim=500000
!===============================
 Integer::I,J
 
 Integer::NC !!Number of Cells of mesh
 Integer::NP !  Number of Existing Points
 Integer::NF !!Number of Faces Constructing Mesh
 Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NFW1,NFW2 !Index of 1st and last Faces on Wall Boundary 
 Integer::NFF1,NFF2 !Index of 1st and Last Faces on Far-Field Boundary
 Integer::NFI1,NFI2 !Index of 1st and Last Faces on Inflow Boundary
 Integer::NFS1,NFS2 !Index of 1st ans Last Faces on Symmetry Boundary
 Integer::NFO1,NFO2 !Index of 1st and Last Faces on Inflow Boundary 
 Integer::NFIF1,NFIF2 !Index of 1st and last Faces on InterFsce Boundary   
 Integer::NS !Number of Steps in rung-kutta method
 Integer::cell !  COUNTER for numBER OF CELLS
 Integer::NBP !Number of Boundary Points
 Integer::Unsteady_Moving !  With this variable first of simulation steady state solution is computed then unsteady moving mesh is computed
 Integer::Test_Case !Number of test case
 Integer::N_Positions  !  Number of Positions to Write Results
 Integer::Position_Num !  Counter of Positions for writing results
 Integer::NSBP         !  Number of Selected Boundary Points
 Integer::Istp         !  index of moving step
 Integer::Stdy_cyc     !  Number Of Steady Cycles
 Integer::UnStdy_cyc   !  Unsteady Cycle Steps
 Integer::II           !!If this term be equal to one steady cycle is completed and unsteady tend to  be started
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::Co    !!  Coefficients of residuals in Rung-kutta steps
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0,V0,W0 !Component of Infinite Flow Velocity
 Real(8)::E0   !Internal Energy of Infinite Flow
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::TT !Total Temprature
 Real(8)::B0 !Sauterland Constant
 Real(8)::RM !Residual of Mass Equation
 Real(8)::U,V,W ! Velocity in X AND Y AND Direction 
 Real(8)::RKco !!  Rung-kutta Coefficients
 Real(8)::Temp   !  value of denominator in the friction of caculating Conservative Values at (N+1)th Time Step by dual time step method
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::DT_Real !Minimum Time Step*Time_Coe
 Real(8)::Real_Time !  Real Time Passed For converging
 Integer::NRKS !Number of Runge Kutta Stages
 Integer::NWrite   !Number of Cycle to Write Results
 Integer::NR   !Number Of Regions of mesh
 Integer::Init !Initialize from Available Data in the file(1) or Infinite Flows(0)
 Real(8)::Total_Time  !  Total time Considerd For Solution By User
 Real(8)::CL !Lift Coefficient
 Real(8)::CD !DRAG Coefficient
 Real(8)::phiWrite !  Local Varible for Counter of Positions for writing results
 Real(8)::Delphi  !!  Differentiation In Phi Direction
 Real(8)::phi      !!  Momentary Angle Of Phase
 Real(8)::phi1     !  Represent Residual of Division phi per 360
 Real(8)::Omega   !Frequency of Flapping Airfoil 
 Real(8)::PI  !Mathematical constant that equals to 3.14159
 Real(8)::X1  !First Point Coordinate of First Line
 Real(8)::H   !  Local Variable that Constitute Horizental Axis of Lift Coefficient Graph
 Real(8)::StdyERmx !Value Of ERror considered for mass equation by user
 Real(8)::Time_Coe  !Time Coefficient for dual time stepping
 Real(8)::UnStdyResMas !Residual of Real Time Iteration
 Real(8)::StdyResMas !  Residual of Imaginary Iteration
 Real(8)::AA  !  Local Variable for Area of each cell
 Real(8)::PlanForm_Area
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Integer,Dimension(1:Dim)::NFace_Cell !Number of Faces Forming a Cell
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:6,1:Dim)::IFace_Cell !Index of Faces Forming a Cell
 Integer,Dimension(1:8,1:Dim)::Corn !Corners point index of Each Cell 
 Real(8),Dimension(1:40)::phi_Write ! !  Phase of Positions to Write Results
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xc,Yc,Zc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::Vol !Volume of each Cell
 Real(8),Dimension(1:Dim)::NX,NY,NZ !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::DT !Time step
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:Dim)::Xo,Yo,Zo !Coordinate of the rotation axis
 Real(8),Dimension(1:Dim)::DelX,DelY,DelZ!!Displacement  in x and y and z Direction
 Real(8),Dimension(1:Dim)::GF !Grid Flux
 Real(8),Dimension(1:3,1:Dim)::FACE_VELOCITY !Velocity of Faces
 Real(8),Dimension(1:5,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WN !Conservative Values at (N)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WNM1 !Conservative Values at (N-1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(1:5,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:5,1:Dim)::Res!!Residual of Equations
 Real(8),Dimension(1:6,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Real(8),Dimension(1:Dim)::SumCP=0.0  !Summation of Averegaed pressure coefficient
 Integer::Naverage=0   !Number of Averaging Steps
!***************************************** Main ********************************************
 PI=4.*Atan(1.)
!Part 1:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3: 
 Call Read_Write_Positions(Test_Case,N_Positions,phi_Write)
 
!Part 4:
 Call MeshBC3D(Dim,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 5:
 Call BoundPointLabeling_3D(Dim,IDS,NR,NFR,BC,FaceType,NBP,IBP)

!Part 6:
 Call FaceOfCell(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 Do Cell=1,NC
    Call PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 End Do
 
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
    
    WN(1,J) = WNP1(1,J)
    WN(2,J) = WNP1(2,J)
    WN(3,J) = WNP1(3,J)
    WN(4,J) = WNP1(4,J) 
    WN(5,J) = WNP1(5,J)
 End Do
 
 !Part 10:
 GF(:) = 0.0
 Face_Velocity(:,:) = 0.0
 
 UnStdyResMas = 100.0  !  Rm = 10.0
 
 UnStdy_cyc = 0   ! Ncyc = 0
 Real_Time = 0.0
 TotTime = 200.0
 
 Unsteady_Moving = 0 ! Unsteady = 0
 
 II=1
 Position_Num = 1
 Istp=1
 
!Part 11:
 Call BC_Wall3D_ALE(Dim,NFW1,NFW2,IDS,GM,P,Face_Velocity,WB)
 Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 
!Part 12:
 Open(100,file='cd.plt')
 Write(100,*)'VARIABLES="t","Cd"'
 
 Open(10,file='CL_Alpha.plt')
 Write(10,*)'VARIABLES="Alpha","CL"'
 
!Part 13:
 Do While(Real_Time < TotTime)
     
   !Part 14:
    UnStdy_cyc = UnStdy_cyc + 1  ! Ncyc=Ncyc+1
     
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
    End Do
    
   !Part 16:
    Call TimSTP_Inviscid3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,DA,CFLx,GM,P,WNP1,WB,DT)
    
   !Part 17:
    DT_Real = 1000.0   !Real_DT
    Do J=1,NC
        If( DT_Real > DT(J) )     DT_Real = DT(J)
    End Do
    DT_Real = DT_Real * Time_Coe
    
   !Part 18:
    Real_Time = Real_Time+DT_Real
    
   !Part 19:
    If(Unsteady_Moving==1) Then!UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

       !Part 20:
        Call Defin_BoundPoint_Displac_3D(Dim,NFW1,NFW2,FaceType,IDS,X,Y,Z,Real_Time,DT_Real,Test_Case,Minf,ALF,Omega,Delx,Dely,Delz) 
        
	   !Part 21:
        Call RBF_GreedyMovingMesh3D(Istp,Dim,NBP,NP,IBP,X,Y,Z,DelX,DelY,DelZ,NSBP,ACTV)
        Istp=Istp+1
        
       !Part 22:
	    Do J=1,NP
            Xo(J)=X(J)
            Yo(J)=Y(J)
            Zo(J)=Z(J)
            X(J)=X(J)+DelX(J)
            Y(J)=Y(J)+DelY(J)
            Z(J)=Z(J)+Delz(J)
        End Do
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        
	   !Part 23:
        Call GridFlux_3D(Dim,NF,IDS,X,Xo,Y,Yo,Z,Zo,FaceType,DT_Real,GF,Face_Velocity)
          
       !Part 24:
        Call GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)
        
       !Part 25:
        Do J=1,NC
            If(Vol(J)<0.) Then
                Print*,'Unable to Move the Boundary'
                Goto 100
            End If
        End Do
        
    End If   !Unsteady_Moving==1  !UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    
   !Part 26:
    StdyResMas = 100.0
    Stdy_cyc = 0
    
   !Part 27:
    Do While(StdyResMas > StdyERmx)
        
       !Part 28:
        Stdy_cyc = Stdy_cyc + 1
	    print*,'Stdy_cyc: ',Stdy_cyc
        
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
            Call ConMeanFlow_ScalarDiss_ALE_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,GF,Con)
           !Call ConMeanFlow_AUSM3D_ALE(Dim,NC,NF1,NF2,NF,GM,IDS,Nx,Ny,Nz,DA,WNP1,WB,P,GF,Con) 
            
           !Part 33:
            Do J=1,NC
                AA = Vol(J)
                Res(1,J) = -Con(1,J)/AA
                Res(2,J) = -Con(2,J)/AA
                Res(3,J) = -Con(3,J)/AA
                Res(4,J) = -Con(4,J)/AA
                Res(5,J) = -Con(5,J)/AA
                
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
                
            End Do !J
            
           !Part 34: 
            Call BC_Wall3D_ALE(Dim,NFW1,NFW2,IDS,GM,P,Face_Velocity,WB)
            Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
            Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
            Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)
            Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
            
        End Do !Ns
        
       !Part 35:
 	    Call ResMass3D(Dim,NC,WNP1,WC,DT,Stdy_cyc,StdyResMas)
        !Print*,'Stdy_cyc: ' , Stdy_cyc , '  StdyResMas:' , StdyResMas
        
    End Do !Steady
    
    !Part 36:
     If( Mod(UnStdy_cyc,NWrite)==0 )Then
        !Part 37:
         Do J=1,NC
             DT(J) = DT_Real
	     End Do
         Call ResMass3D(Dim,NC,WNP1,WN,DT,UnStdy_cyc,UnStdyResMas)
         Print*,UnStdy_cyc,UnStdyResMas,Real_Time
         
         !Part 38:
         If( UnStdyResMas < -3.0 .AND. II==1 ) Then
             II=2
             Unsteady_Moving=1
             UnStdy_cyc = 0
             Real_Time = 0.0
         End If
         
         !Part 39:
         Print*, ' Alfa = ' , ALF
         
         !Part 40:
         If (Test_Case <= 4) Then
             X1 = ALF
         Else If(Test_Case == 5) Then
             H = 0.0125*(sin(0.8*Real_Time))
             X1 = H
         Else If(Test_Case == 6) Then
             X1 = omega*Real_Time
         Else If(Test_Case == 7) Then
             X1 = Real_Time * Minf
         Else If(Test_Case == 8) Then
             X1 = Real_Time
         End If
         
         Call PressLiftDragCo_3D(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,NZ,DA,Unsteady_Moving,PlanForm_Area,CL,CD)
         Print*,PlanForm_Area
         Write(10,*)  X1,CL
         Write(100,*) X1,CD
         
     End If
     
     !Part 41:
     If(Unsteady_Moving==1) Then!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

         phiWrite = phi_Write(Position_Num)
         phi  = omega*Real_Time*(180.0/PI)
         phi1 = Mod(phi,360.0)
         Delphi = 1.2*omega*DT_Real*(180.0/PI)
         If( phi1 > (phiWrite-Delphi) .AND. phi1 < (phiWrite+Delphi) ) Then
             Position_Num = Position_Num + 1
             
             Call Write_VelocityContour3D(Dim,NC,NP,X,Y,Z,Corn,WNP1)
             Call Write_CP3DUnsteady(Dim,NFW1,NFW2,IDS,FaceType,X,Y,Z,Minf,P,GM,Naverage,SumCP)
             
             If( N_Positions == Position_Num-1 ) Position_Num = 1
         End If

     End If !Unsteady_Moving==1!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         
 End Do !Do While
!*********************************************************************************************
 Close(100)
 Close(10)
!Part 42:
100 End 
!###########################################################################################


