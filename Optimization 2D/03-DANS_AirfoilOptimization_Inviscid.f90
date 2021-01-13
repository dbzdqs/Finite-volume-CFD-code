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
!// This code are designed to optimize airfoil shape which the flow solver of this code is an inviscid //!
!// flow solver. The moving mesh method for updating mesh as the geometry modified through optimization//!
!// algorithm is RBF moving mesh method.                                                               //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_AirfoilOptimization_Inviscid
 Implicit None
!********************************************************************************************* 
!===============================
 Integer,Parameter::Dim = 10000
 Integer,Parameter::SwarmSize = 5		! Number of Swarms considered 
 Integer,Parameter::Iterations =300	    ! Particle swarm optimization Maximum iterations number
 Real(8),Parameter::ConvergERROR = 1e-2	! Particle swarm optimization Convergance condition
!===============================

 Integer::I,J,S,K,P1
 Integer::NRKS !Number of Runge Kutta Stages
 Integer::NWrite   !Number of Cycle to Write Results
 Integer::Init !Initialize from Available Data in the file(1) or Infinite Flows(0)
 Integer::NP  !Number of Existing Points
 Integer::NC !Number of Cells of mesh
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NR   !Number Of Regions of mesh
 Integer::NBP !Number of Boundary Points
 Integer::UpRegion
 Integer::LwRegion
 Integer::DimU
 Integer::DimL
 Integer::BSOrder_Up
 Integer::BSOrder_Lw
 Integer::NPtCurv_Up,NPtCurv_Lw
 Integer::Check
 Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NFW1,NFW2 !Index of 1st and last Faces on Wall Boundary 
 Integer::NFF1,NFF2 !Index of 1st and Last Faces on Far-Field Boundary
 Integer::NFI1,NFI2 !Index of 1st and Last Faces on Inflow Boundary 
 Integer::NFS1,NFS2 !Index of 1st ans Last Faces on Symmetry Boundary
 Integer::NFO1,NFO2 !Index of 1st and Last Faces on Inflow Boundary
 Integer::NFIF1,NFIF2 !Index of 1st and last Faces on InterFsce Boundary
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points													  
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Allocatable,Dimension(:)::IBP_Lw,IBP_Up
 Real(8)::Zita_TE
 Real(8)::TE_Thick
 Real(8)::Shap
 Real(8)::Xo,Yo
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0   !Infinite Flow Velocity in X Direction
 Real(8)::V0   !Infinite Flow Velocity in Y Direction
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::U,V 
 Real(8)::CL !Lift Coefficient
 Real(8)::CD !Drag Coefficient
 Real(8)::ObjGbest
 Real(8)::ObjectiveFun
 Real(8) ::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TT !Total Temprature
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient
 Real(8),Dimension(1:Dim)::Xc,Yc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::A !Area of each cell
 Real(8),Dimension(1:Dim)::NX,NY !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::DelX,DelY
 Real(8),Dimension(1:4,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,Dimension(1:4,1:Dim)::Corn !Corners point index of Each Cell 
 Real(8),Allocatable,Dimension(:)::Pol_Coeff_Up
 Real(8),Allocatable,Dimension(:)::Pol_Coeff_Lw
 Real(8),Allocatable,Dimension(:)::X_Up
 Real(8),Allocatable,Dimension(:)::Y_Up
 Real(8),Allocatable,Dimension(:)::X_Lw
 Real(8),Allocatable,Dimension(:)::Y_Lw
 Real(8),Allocatable,Dimension(:)::S_Lw
 Real(8),Allocatable,Dimension(:)::S_Up
 Real(8),Allocatable,Dimension(:)::Ynew_Up
 Real(8),Allocatable,Dimension(:)::Ynew_Lw
 Real(8),Allocatable,Dimension(:)::ShapeFunc_Coe
 Real(8),Allocatable,Dimension(:)::gBest
 Real(8),Allocatable,Dimension(:)::GeoCoeff
 Real(8),Allocatable,Dimension(:)::Pmax
 Real(8),Allocatable,Dimension(:)::Pmin
 Real(8),Allocatable,Dimension(:)::ObjPbest
 Real(8),Allocatable,Dimension(:)::ObjFuncS
 Real(8),Allocatable,Dimension(:,:)::BMI_Up
 Real(8),Allocatable,Dimension(:,:)::BMI_Lw
 Real(8),Allocatable,Dimension(:,:)::Pos
 Real(8),Allocatable,Dimension(:,:)::Velocity
 Real(8),Allocatable,Dimension(:,:)::pBest
 Real(8),Allocatable,Dimension(:,:)::CheckMatrix  
 Logical:: Converged = .False.
!**************************************Main*****************************************************************	
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 
!Part 2:
 Call Read_CST_Parameters(BSOrder_Up,BSOrder_Lw,UpRegion,LwRegion)
!============================================ Airfoli Parametization ================================
!Part 3:
 DimU=NFR(UpRegion)+1
 DimL=NFR(LwRegion)+1

!Part 4:
 Allocate( Pol_Coeff_Lw(0:BSOrder_Lw)            , Pol_Coeff_Up(0:BSOrder_Up),&          
           BMI_Lw(1:BSOrder_Lw+1,1:BSOrder_Lw+1) , BMI_Up(1:BSOrder_Up+1,1:BSOrder_Up+1),& 
		   IBP_Lw(1:DimL)                        , IBP_Up(1:DimU),&
		   X_Lw(1:DimL)                          , X_Up(1:DimU),&
		   Y_Lw(1:DimL)                          , Y_Up(1:DimU),&
		   S_Lw(1:DimL)                          , S_Up(1:DimU),&
		   Ynew_Lw(1:DimL)                       , Ynew_Up(1:DimU),&
           ShapeFunc_Coe(1:BSOrder_Up+1+BSOrder_Lw+1) )


!Part 5:
 Call CST_PreProcc(Dim,DimU,DimL,NP,NR,NFR,IDS,UpRegion,LwRegion,X,Y,&
                   NPtCurv_Lw,NPtCurv_Up,Zita_TE,TE_Thick,IBP_Lw,IBP_Up,X_Up,Y_Up,X_Lw,Y_Lw)

!Part 6:
 Call Shape_Function(DimU,DimL,NPtCurv_Up,X_Up,Y_Up,NPtCurv_Lw,X_Lw,Y_Lw,Zita_TE,S_Up,S_Lw)

!Part 7:
 Call Bernstein_Matrix_Inverse(BSOrder_Up+1,BMI_Up)
 Call Bernstein_Matrix_Inverse(BSOrder_Lw+1,BMI_Lw)

!Part 8:
 Call Poly_Fit(NPtCurv_Up,X_Up,S_Up,BSOrder_Up,Pol_Coeff_Up)
 Call Poly_Fit(NPtCurv_Lw,X_Lw,S_Lw,BSOrder_Lw,Pol_Coeff_Lw)

!Part 9:
 Do I=1,BSOrder_Up+1

    Shap=0.0
    Do J = 1, BSOrder_Up+1
       Shap = Shap + BMI_Up(I,J)*Pol_Coeff_up(J-1)				
    End Do

  ShapeFunc_Coe(I)= Shap
 End Do

!Part 10:
 Do I=1,BSOrder_Lw+1

    Shap=0
    Do J = 1, BSOrder_Lw+1
       Shap = Shap + BMI_Lw(I,J)*Pol_Coeff_Lw(J-1)
    End Do

    ShapeFunc_Coe(I+BSOrder_Up+1)= Shap
 End Do

 Deallocate( Pol_Coeff_Lw , Pol_Coeff_Up , BMI_Lw , BMI_Up )
!========================================= End of Airfoli Parametization ================================

!========================================= Flow Solver Initialization ===================================
!Part 11:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)

!Part 12:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 13:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 14:
 Call InitMeanFlow_Inviscid(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,WNP1)

!Part 15:
 Do J=1,NC
    U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do
!====================================== End of Flow Solver Initialization ================================

!Part 16:
Call BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)

!Part 17:
Allocate( Pos(1:SwarmSize,1:BSOrder_Up+1+BSOrder_Lw+1),Velocity(1:SwarmSize,1:BSOrder_Up+1+BSOrder_Lw+1),&
		   pBest(1:SwarmSize,1:BSOrder_Up+1+BSOrder_Lw+1),gBest(1:BSOrder_Up+1+BSOrder_Lw+1),&
		   Pmax(1:BSOrder_Up+1+BSOrder_Lw+1),Pmin(1:BSOrder_Up+1+BSOrder_Lw+1),ObjpBest(1:SwarmSize),&
		   CheckMatrix(1:SwarmSize,1:BSOrder_Up+1+BSOrder_Lw+1),ObjFuncS(1:SwarmSize),GeoCoeff(1:BSOrder_Up+1+BSOrder_Lw+1) )

 Open(15,File='Check.dat')
 Open(16,File='OptimizedShape.dat')
 Open(55,File='Converge Check.dat')

!Part 18:  
 Call PSO_Initial(SwarmSize,BSOrder_Up+1+BSOrder_Lw+1,ShapeFunc_Coe,Pmax,Pmin,Pos,Velocity,pBest)

!Part 19:
 Do I=1,Swarmsize
   !minmizing
    ObjPbest(I) = 10000000
   !Maximizing 
   !ObjPbest(I) = 0.0000001
 End Do

!minimizing
 ObjGbest = 10000000
!Maximizing 
!ObjGbest = 0.00000001

!Part 20:
 K = 0				
 Do While ( (K <= Iterations).And.(.Not.Converged) )

    K = K+1
    
    Print*,'Optimization Iterations:',K 

   !Part 21:
	Do S = 1 ,SwarmSize

	    Write(15,'(a,I4,I4)') 'K - S - GeoCoeff', K, S
   	    Do I = 1,BSOrder_Up+1+BSOrder_Lw+1
	       GeoCoeff(I) = 	 Pos(S,I)
	       Write(15,'(F10.6)')  GeoCoeff(I)
        End Do

       !Part 22:
    	Call CST_InversToShap(NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,GeoCoeff,Zita_te,X_Up,X_Lw,Ynew_Up,Ynew_Lw)
        
       !Part 23:
        Call Update_Mesh(Dim,DimU,DimL,NPtCurv_Up,IBP_UP,Ynew_Up,Y_Up,NPtCurv_Lw,IBP_Lw,Ynew_Lw,Y_Lw,NBP,NP,NC,NF,&
                        NF1,NF2,IBP,X,Y,IDS,DelX,DelY,Xc,Yc,NX,NY,DA,A)
     
       !Part 24:
        Call Solver_DANSInviscid(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,Minf,ALF,ERmx,CFLx,NRKS,NWrite,RKJ,&
							        NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
							        Xc,Yc,NX,NY,DA,A,GM,R0,P0,C0,U0,V0,&
							        WNP1,P,WB) 
    
	   !Part 25: 
        Call PressLiftDragCo(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,CL,CD)

        ObjectiveFun = abs(Cd/Cl)

	    Write(15,'(a,I4)') 'CD  CL',S 							
        Write(15,'(F10.6,F10.6)') CD,CL 

        Write(15,'(a,I4,I4)') 'K - S - ObjectiveFun' , K , S
        Write(15,'(F10.6)') ObjectiveFun

        Print*,' Swarm :' , S ,'ObjectiveFun' , ObjectiveFun 	 

	    ObjFuncS(S) =  ObjectiveFun

	   !Part 26:
       !For Maximizing >
       !For minimizing <    
	    If( ObjectiveFun < ObjPbest(S) )Then
         Do I = 1,BSOrder_Up+1+BSOrder_Lw+1
	        pBest(S,I) = Pos(S,I)
         End Do
         ObjPbest(S) = ObjectiveFun
        End If
	
        Write(15,'(a,I4)') 'S - ObjPbest' , S
        Write(15,'(F10.6)') ObjPbest(S)

       !Part 27:
       !For Maximizing >
       !For minimizing <
        If(ObjPbest(S) <  ObjGbest )Then
	       Do I = 1,BSOrder_Up+1+BSOrder_Lw+1
              gBest(I) = pBest(S,I)
           End Do
           ObjGbest = ObjPbest(S)
        End If
			   
        Write(15,'(a,I4)') ' k - ObjGbest' , K
        Write(15,'(F10.6)') ObjGbest

       !Part 28:
        Do I = 1, NPtCurv_Up
           Y_Up(I) = Ynew_Up(I)
        End Do
    
	    Do I = 1, NPtCurv_Lw
           Y_Lw(I) = Ynew_Lw(I)
	    End Do  

    End Do

   !Part 29:
    Call PSO_Update(K,Iterations,SwarmSize,BSOrder_Up+1+BSOrder_Lw+1,Pmax,Pmin,pBest,gBest,ObjFuncS,ObjGbest,Pos,Velocity)

   !Part 30:
 	Check = 0
    Do S = 1,SwarmSize
	   If (ABS(ObjFuncS(S) - ObjGbest)<= ConvergERROR) Check = Check + 1
   	End Do
  
    If (Check == SwarmSize ) converged = .True.

   !Part 31:
    Write(55,'(I4,F10.6,F10.6,F10.6,F10.6,F10.6)') K , ObjFuncS(1) , ObjFuncS(2), ObjFuncS(3), ObjFuncS(4), ObjFuncS(5)

 End Do
    
!Part 32:
 Write(15,'(a)') 'gBest'
  
 Do I = 1,BSOrder_Up+1+BSOrder_Lw+1
    GeoCoeff(I) = gBest(I)
    Write(15,'(F10.6)')  GeoCoeff(I)
 End Do

 Call CST_InversToShap(NPtCurv_Up,NPtCurv_Lw,BSOrder_Up,BSOrder_Lw,GeoCoeff,Zita_te,X_Up,X_Lw,Ynew_Up,Ynew_Lw)


 Write(16,'(a)') 'optimized CST_Up'
 Do I= 1, DimU
    Write(15,'(F10.6,F10.6)')  X_Up(I) , Ynew_Up(I)
 End Do

 Write(16,'(a)') 'Optimized CST_Low'
 Do I= 1, DimL
    write(15,'(F10.6,F10.6)')  X_Lw(I) , Ynew_Lw(I)
 End Do

 Call Update_Mesh(Dim,DimU,DimL,NPtCurv_Up,IBP_UP,Ynew_Up,Y_Up,NPtCurv_Lw,IBP_Lw,Ynew_Lw,Y_Lw,NBP,NP,NC,NF,&
                  NF1,NF2,IBP,X,Y,IDS,DelX,DelY,Xc,Yc,NX,NY,DA,A)

 Call Solver_DANSInviscid(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,Minf,ALF,ERmx,CFLx,NRKS,NWrite,RKJ,&
                             NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
							 Xc,Yc,NX,NY,DA,A,GM,R0,P0,C0,U0,V0,&
							 WNP1,P,WB)  

 Call PressLiftDragCo(Dim,NFW1,NFW2,GM,Minf,ALF,WB,NX,NY,DA,CL,CD)

 Write(16,*) '=============final results============='	
 Write(16,'(a,F10.6)') 'final CL' , CL		
 Write(16,'(a,F10.6)') 'final Cd' , Cd
 Write(16,'(a,F10.6)') 'Final ObjGbest' , ObjGbest

 Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
 Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)

!*********************************************************************************************
 End Program DANS_AirfoilOptimization_Inviscid
!###########################################################################################
