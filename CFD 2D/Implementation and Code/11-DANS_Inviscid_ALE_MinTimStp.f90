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
!// 1-Dimension:	                            2D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Invicid/ALE                                            //!
!// 7-Convection Discretization Scheme:	        AUSM                                                   //!
!//                                             AUSM+                                                  //!
!//                                             Ausm+Up                                                //!
!//                                             HLL                                                    //!
!//                                             HLLC                                                   //!
!//                                             MatrixDiss                                             //!
!//                                             Roe_ENT                                                //!
!//                                             Roe_JCB                                                //!
!//                                             RoeEC                                                  //!
!//                                             ScalarDiss                                             //!
!//                                             KEP                                                    //!
!//                                             CUSM Zha                                               //!
!//                                             CUSP98                                                 //!
!// 8-Convection Term Discretization Accuracy:	First order                                            //!
!// 9-Transient Term Discretization Scheme:     Rung-Kuta, Minimum Time Stepping                       //!
!// 10-Steady/Unsteady: 	                    Unsteady                                               //!
!// 11-Gradient Calculation Scheme:	            non                                                    //!
!// 12-Turbulence Model:	                    non                                                    //!
!// 13-Moving Mesh Method:	                    RBF/Linear Spring                                      //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Inviscid_ALE_MinTimStp
 Implicit None
!===============================
 Integer,Parameter::Dim=50000
!===============================
 Integer::I,J
 Integer::NS    !  Number of Steps in rung-kutta method
 Real(8)::U    !  Velocity in x Direction
 Real(8)::V    !  Velocity in Y Direction
 Real(8)::RKco !  Rung-kutta Coefficients
 Real(8)::co   !  Coefficients of residuals in Rung-kutta steps 
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
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0   !Infinite Flow Velocity in X Direction
 Real(8)::V0   !Infinite Flow Velocity in Y Direction
 Real(8)::E0 !Internal Energy of Infinite Flow
 Real(8)::H0 !Enthalpy of Infinite Flow
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::RM !Residual of Mass Equation
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(1:4,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:4,1:Dim)::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(1:4,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xc,Yc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::A !Area of each cell
 Real(8),Dimension(1:Dim)::NX,NY !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::DT !Time step
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:5,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Integer,Dimension(1:4, 1:Dim)::Corn !Corners point index of Each Cell 
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TT !Total Temprature
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Integer::NBP !Number of Boundary Points
 Integer::Unsteady !  With this variable first of simulation steady state solution is computed then unsteady moving mesh is computed
 Integer::Test_Case !!Number of test case
 Integer::N_Positions !Number of Positions to Write Results
 Integer::Position_Num  !Counter of Positions for writing results
 Real(8)::DT_REAL !Minimum Time Step*Time_Coe
 Real(8)::Real_Time !  Real Time Passed For converging
 Real(8)::CL !Lift Coefficient
 Real(8)::CD !Drag Coefficient
 Real(8)::H !!  Local Variable that Constitute Horizental Axis of Lift Coefficient Graph
 Real(8)::phiWrite!!Phase of Positions to Write Results
 Real(8)::Delphi!!  Differentiation In Phi Direction
 Real(8)::phi  !!  Momentary Angle Of Phase
 Real(8)::phi1  !!  Represent Residual of Division phi per 360
 Real(8)::Omega  !Frequency of Flapping Airfoil
 Real(8)::PI     !Mathematical constant that equals to 3.14159
 Real(8)::X1     !First Point Coordinate of First Line !First Point Coordinate of First Line
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Real(8),Dimension(1:Dim)::DelX,DelY !!Displacement  in x and y Direction
 Real(8),Dimension(1:Dim)::GF !Grid Flux
 Real(8),Dimension(1:2,1:Dim)::FACE_VELOCITY !Velocity of Faces
 Real(8),Dimension(1:40)::phi_Write!Phase of Positions to Write Results
!************************************************** Main ************************************************
 PI=4.*Atan(1.)
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 
!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3: 
 Call Read_Write_Positions(Test_Case,N_Positions,phi_Write)
  
!Part 4:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 5:
 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
!Part 6:
 Call PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,IDS,Corn)
 
!Part 7:
 Call BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)
 
!Part 8:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)
 
!Part 9:
 Call InitMeanFlow_Inviscid(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,WNP1)
     
!Part 10: 
 Do J=1,NC
     U = WNP1(2,J)/WNP1(1,J)
     V = WNP1(3,J)/WNP1(1,J)
     P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do
 
!Part 11:
 GF = 0.0
 Face_Velocity = 0.0
 Rm = 10.0
 Ncyc = 0
 Real_Time = 0.0
 Unsteady = 0
 Position_Num = 1
 
!Part 12:
 Call BC_InvisWall_ALE(Dim,NFW1,NFW2,IDS,GM,P,NX,NY,DA,WNP1,Face_Velocity,WB)
!Call BC_Wall_ALE(Dim,NFW1,NFW2,IDS,GM,WNP1,P,Face_Velocity,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
 
!Part 13:
 Open(100,file='cd.plt')
 Write(100,*)'VARIABLES="t","Cd"'
 
 Open(10,file='CL_Alpha.plt')
 Write(10,*)'VARIABLES="Alpha","CL"'
 
!Part 14:
 Do While(Real_Time < TotTime)
     
    !Part 15:
     Ncyc=Ncyc+1
     
    !Part 16:
     Do J=1,NC
         WC(1,J) = WNP1(1,J)
         WC(2,J) = WNP1(2,J)
         WC(3,J) = WNP1(3,J)
         WC(4,J) = WNP1(4,J)   
     End Do
     
    !Part 17:
     Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)
     
    !Part 18:
     If(Unsteady==1) Then
       
        !Part 19:
         DT_Real = 1000.0
         Do J=1,NC
             If( DT_Real > DT(J) )     DT_Real = DT(J)
         End Do
	     DT=DT_Real
         
        !Part 20:
	     Real_Time = Real_Time+DT_Real
        
        !Part 21:
         Call Defin_BoundPoint_Displac(Dim,NFW1,NFW2,IDS,X,Y,Real_Time,DT_Real,Test_Case,Minf,ALF,Omega,Delx,Dely) 
         
        !Part 22:
         Call RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)
         
        !Part 23:
	     Do J=1,NP
             Xo(J)=X(J)
             Yo(J)=Y(J)
             X(J)=X(J)+DelX(J)
             Y(J)=Y(J)+DelY(J)
         End Do
         
        !Part 24:
         Call GridFlux(Dim,NF1,NF2,NFW1,NFW2,IDS,X,Xo,Y,Yo,DT_Real,GF,Face_Velocity)
         
        !Part 25:
         Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A) !,A1
         
        !Part 26:
         Do j=1,NC
             If(A(J)<0.) Then
                 Print*,'negative Volume Unable to Move the Boundary'
                 Goto 100
             End If
         End Do
       
     End If    
     
    !Part 27:
     Do NS=1,NRKS
         
        !Part 28:
         RKco=1.0/(NRKS-NS+1)
         
        !Part 29:
         Call ConMeanFlow_ScalarDiss_ALE(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,K2,K4,WNP1,WB,P,GF,Con)
        !Call ConMeanFlow_AUSM_ALE(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,GF,WNP1,WB,P,Con)
         
        !Part 30:
         Do J=1,NC
             Co = RKco*DT(J)/A(J)
             WNP1(1,J) = WC(1,J) - Co* Con(1,J) 
             WNP1(2,J) = WC(2,J) - Co* Con(2,J)
             WNP1(3,J) = WC(3,J) - Co* Con(3,J)
             WNP1(4,J) = WC(4,J) - Co* Con(4,J)
             
             U = WNP1(2,J) / WNP1(1,J)
             V = WNP1(3,J) / WNP1(1,J)
             P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
         End Do
         
       !Part 31:
        Call BC_InvisWall_ALE(Dim,NFW1,NFW2,IDS,GM,P,NX,NY,DA,WNP1,Face_Velocity,WB)
        ! Call BC_Wall_ALE(Dim,NFW1,NFW2,IDS,GM,WNP1,P,Face_Velocity,WB)
         Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
         Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
         Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
         Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
         
     End Do !Ns	
     
    !Part 32:
     If( Mod(Ncyc,NWrite)==0 ) Then
         Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
         Print*,Ncyc,Rm,Real_Time
         
         If(Rm < ERmx) Unsteady=1
         
         Print*, ' Alfa = ' , ALF
         
         !Part 33:
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
         
         Call PressLiftDragCo(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,CL,CD)
         Write(10,*)  X1,CL
         Write(100,*) X1,CD
     End If
     
    !Part 34:
     If(Unsteady==1) Then

         phiWrite = phi_Write(Position_Num)
         phi  = omega*Real_Time*(180.0/PI)
         phi1 = Mod(phi,360.0)
         Delphi = 1.2*omega*DT_Real*(180.0/PI)
         If( phi1 > (phiWrite-Delphi) .AND. phi1 < (phiWrite+Delphi) ) Then
             Position_Num = Position_Num + 1
              Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
              Call Write_ConservativeVariables(Dim,NC,WNP1)
	          Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
             If( N_Positions == Position_Num-1 ) Position_Num = 1
         End If

     End If !Unsteady==1
     
 End Do !Do While
 
 Close(100)
 Close(10)

 100 End 
!########################################################################################################