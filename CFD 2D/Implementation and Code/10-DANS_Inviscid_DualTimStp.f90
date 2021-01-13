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
!// This Code designed for solving unsteady and inviscid form of two dimensional Euler equations. One  //!
!// of the main advantages of this code is that it enjoy several subroutines that assists user to      //!
!// follow code step by step. Mesh and Settings have external files that read through their            //!
!// subroutines in the first step of code. In addition dual time stepping scheme used for time         //!
!// marching. Main features of this code is as following:                                              //!
!// 1-Dimension:	                            2D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Invicid                                                //!
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
!// 9-Transient Term Discretization Scheme:     Dual Time Stepping                                     //!
!// 10-Steady/Unsteady: 	                    Unsteady                                               //!
!// 11-Gradient Calculation Scheme:	            non                                                    //!
!// 12-Turbulence Model:	                    non                                                    //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Inviscid_DualTimStp
 Implicit None
!===============================
 Integer,Parameter::Dim=120000
!===============================

 Integer::I,J
 Integer::NS    !  Number of Steps in rung-kutta method
 Real(8)::U    !  Velocity in x Direction
 Real(8)::V    !  Velocity in Y Direction
 Real(8)::RKco !  Rung-kutta Coefficients
 Real(8)::co   !  Coefficients of residuals in Rung-kutta steps
 Real(8)::Temp  !  value of denominator in the friction of caculating Conservative Values at (N+1)th Time Step by dual time step method
 Real(8)::AA    !  Local Variable for Area of each cell
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
 Real(8),Dimension(1:4,1:Dim)::WN !Conservative Values at (N)th Time Step
 Real(8),Dimension(1:4,1:Dim)::WNM1 !Conservative Values at (N-1)th Time Step
 Real(8),Dimension(1:4,1:Dim)::Res!!Residual of Equations
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
 Real(8)::Real_Time !  Real Time Passed For converging
 Real(8)::DT_REAL !Minimum Time Step*Time_Coe
 Integer::Stdy_cyc !  Number Of Steady Cycles
!************************************************** Main **************************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)
 
!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 4:
 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
!Part 5:
 Call PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,IDS,Corn)
 
!Part 6:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 7:
 Call InitMeanFlow_Inviscid(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,WNP1)
     
!Part 8: 
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

!Part 9:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 10:
 Ncyc = 0
 Rm   = 10.0
 Real_Time = 0.0

!Part 11:
 Do While(Real_Time < TotTime)

   !Part 12:
    Ncyc=Ncyc+1
     
   !Part 13:
    Do J=1,NC
         WNM1(1,J) = WN(1,J)
         WNM1(2,J) = WN(2,J)
         WNM1(3,J) = WN(3,J)
         WNM1(4,J) = WN(4,J) 
         
         WN(1,J) = WNP1(1,J)
         WN(2,J) = WNP1(2,J)
         WN(3,J) = WNP1(3,J)
         WN(4,J) = WNP1(4,J)  
    End Do

   !Part 14:
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)
     
    !Part 15:
     DT_Real = 1000.0
     Do J=1,NC
         IF( DT_Real > DT(J) ) DT_Real = DT(J)
     End Do
     DT_Real = DT_Real * Time_Coe
     
    !Part 16:
     Real_Time = Real_Time+DT_Real

    !Part 17:
     Rm = 100.0
     Stdy_cyc = 0
     
    !Part 18:
     Do While(Rm > ERmx .or. Stdy_cyc<2000)
     
        !Part 19:
         Stdy_cyc = Stdy_cyc + 1
	     
        !Part 20:
         Do J=1,NC
             WC(1,J) = WNP1(1,J)
             WC(2,J) = WNP1(2,J)
             WC(3,J) = WNP1(3,J)
             WC(4,J) = WNP1(4,J)   
         End Do
         
        !Part 21:
         Do NS=1,NRKS
            
            !Part 22:
	         RKco=RKJ(NS)
             
            !Part 23:
             Call  ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
             
            !Part 24:
             Do J=1,NC
                 AA = A(J)
                 Res(1,J) = -Con(1,J)/AA
                 Res(2,J) = -Con(2,J)/AA
                 Res(3,J) = -Con(3,J)/AA
                 Res(4,J) = -Con(4,J)/AA
                 
                 Temp = 2*DT_Real+3*RKco*DT(J)
                 WNP1(1,J) = ( 2*WC(1,J)*DT_Real + RKco*DT(J)*(4*WN(1,J) - WNM1(1,J) + 2*Res(1,J)*DT_Real) ) / Temp
                 WNP1(2,J) = ( 2*WC(2,J)*DT_Real + RKco*DT(J)*(4*WN(2,J) - WNM1(2,J) + 2*Res(2,J)*DT_Real) ) / Temp
                 WNP1(3,J) = ( 2*WC(3,J)*DT_Real + RKco*DT(J)*(4*WN(3,J) - WNM1(3,J) + 2*Res(3,J)*DT_Real) ) / Temp
                 WNP1(4,J) = ( 2*WC(4,J)*DT_Real + RKco*DT(J)*(4*WN(4,J) - WNM1(4,J) + 2*Res(4,J)*DT_Real) ) / Temp	
                 
                 U = WNP1(2,J)/WNP1(1,J)
                 V = WNP1(3,J)/WNP1(1,J)
                 P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
                 
             End Do !J
             
            !Part 25:
             Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
             Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
             Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
             Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
             Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)
             
         End Do !Ns	
         
	    !Part 26:
 	     Call ResMass(Dim,NC,WNP1,WC,DT,Stdy_cyc,Rm)
         Print*,'Stdy_cyc: ' , Stdy_cyc , '  ResMas:' , Rm
         
     End Do !Steady
     
	!Part 27:
     Print*,Ncyc,Real_Time
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
      Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
      Call Write_ConservativeVariables(Dim,NC,WNP1)
	  Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
     End If
     
 End Do

!**********************************************************************************************************
 End 
!##########################################################################################################