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
!// This Code designed for Analyzing Turbulence flow by solving two dimensional Navier-Stocks          //!
!// equations. One of the overt advantage of this code is that it enjoy several subroutines that       //!
!// assists user to follow code step by step. Mesh and Settings have external files that read through  //!
!// their subroutines in the first step of code. The code uses Spalart-Allmaras turbulent model for    //! 
!// calculating turbulence viscosity Main features of this code is as following:                       //!
!// 1-Dimension:	                            2D                                                     //!
!// 2-Type Of Mesh:	                            Unstructured                                           //!
!// 3-Data Structure of Mesh:                   Edge Base                                              //!
!// 4-Data Structure of Solver:	                Cell Centered                                          //!
!// 5-Flow Solver algorithm:	                Density Base                                           //!
!// 6-Flow Regime:	                            Turbulent                                              //!
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
!// 9-Transient Term Discretization Scheme:     Explicit(Rung-Kutta)                                   //!
!// 10-Steady/Unsteady: 	                    Steady and Unsteady                                    //!
!// 11-Gradient Calculation Scheme:	            Green-Gauss                                            //!
!// 12-Turbulence Model:	                    Spalart-Allmaras                                       //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_1eqTurbModel2D
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
 Real(8)::Temp  !  Temprature
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
 Real(8)::T0   !Temperature of Infinite Flow
 Real(8)::E0 !Internal Energy of Infinite Flow
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
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
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
 Integer,Dimension(1:4, 1:Dim)::Corn !Corners point index of Each Cell 
 Real(8),Dimension(1:Dim)::TurbQ !Turbulence Quantity (i.e. k in turbulence model)
 Real(8)::Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1 !Spalart-Allmaras Parameters
 Real(8)::Nuhat0  !Turbulence Variable in Spalart-Allmaras
 Real(8),Dimension(1:Dim)::WTNP1 !Turbulence Variables at new time step
!***************************************** Main ********************************************
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
 Call WallDist(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,INW,DW)
 
!Part 7:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 8:
 Call InitMeanFlow(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,T0,Mu0,B0,WNP1)
 
!Part 9: 
 Do J=1,NC     
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
    
   !Part 10:
    Temp = GM*P(J)/WNP1(1,J) 
    
   !Part 11:
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)  
    
 End Do
 
!Part 12:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 13:
 Call SA_Init(Dim,NC,Mu0,R0,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mut0,Nuhat0,WTNP1,Mut) 
 
!Part 14:
 Ncyc = 0
 Rm   = 10.0

!Part 15:
 Do While(Rm > ERmx)

   !Part 16:
    Ncyc=Ncyc+1
    
   !Part 17:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)   
    End Do

   !Part 18:
	Call TimSTP_Turb(Dim,NC,NF,NF1,NF2,IDS,NX,NY,DA,A,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)
    
    !Part 19:
    Do NS=1,NRKS
	   RKco=RKJ(NS)
       
      !Part 20:
	   Call ConMeanFlow_AUSM_PlusUP(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Minf,WNP1,WB,P,Con)

      !Part 21:
       Call VelTemp_GradFace(Dim,NC,NF1,NF2,NFW1,NFW2,NF,NP,IDS,X,Y,Xc,Yc,WNP1,WB,GM,P,DUX,DUY,DVX,DVY,DTX,DTY)

      !Part 22:
       TurbQ(:)=0.0
	   Call DifMeanFlow_TurbNoWallFu(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,GM,PrL,PrT,NX,NY,MR,Mu,Mut,&
                                     WNP1,TurbQ,WB,DUX,DUY,DVX,DVY,DTX,DTY,Dif)
      !Part 23:
       Do J=1,NC
		  Co = RKco*DT(J)/A(J)
          WNP1(1,J) = WC(1,J) - Co*( Con(1,J)          )
		  WNP1(2,J) = WC(2,J) - Co*( Con(2,J)+Dif(2,J) )
          WNP1(3,J) = WC(3,J) - Co*( Con(3,J)+Dif(3,J) )
          WNP1(4,J) = WC(4,J) - Co*( Con(4,J)+Dif(4,J) )

         !Part 24:
	      U    = WNP1(2,J)/WNP1(1,J)
          V    = WNP1(3,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V)) !

         !Part 25:
          Temp = GM*P(J)/WNP1(1,J) 
          
         !Part 26:
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         

       End Do

      !Part 27: 
       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
       Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

     End Do !Ns	

    !Part 28:
     Call SA_Main(Dim,NFW1,NFW2,NF1,NF2,NFF1,NFF2,NFS1,NFS2,NF,NC,NP,IDS,NX,NY,X,Y,Xc,&
                  Yc,A,Dw,NRKS,MR,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mu0,Nuhat0,Mu,DT,WB,WNP1,WTNP1,Mut)
     
    !Part 29:
 	 Call ResMass(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,Ncyc,Rm

    !Part 30:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm
      Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
      Call Write_CF(Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu)
      Call Write_ConservativeVariables(Dim,NC,WNP1)
	  Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mu,"MolVis.plt")
      Call Write_ScalarContour2D(Dim,NC,NP,X,Y,Corn,Mut,"TubVis.plt")
      
     !just for Flat Plate 
      !Call Write_FlatPlateAnalyticData2D(Dim,NFW1,NFW2,NC,IDS,X,Xc,Yc,Minf,Rinf,WNP1,Mu,DUY)
	 End If

 End Do !Do While
!*********************************************************************************************
 End 
!###########################################################################################


