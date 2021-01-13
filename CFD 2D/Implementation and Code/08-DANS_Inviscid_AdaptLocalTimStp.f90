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
!// This Code designed for solving Inviscid form of two dimensional Navier-Stocks equations that known //!
!// as Euler equations. One of the overt advantage of this code is that it enjoy several subroutines   //!
!// that assists user to follow code step by step. Mesh and Settings have external files that read     //!
!// through their subroutines in the first step of code. The main objective of this code is adaptive   //!
!// local time stepping scheme used for time marching. Main features of this code is as following:     //!                                                                                    //!
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
!// 9-Transient Term Discretization Scheme:     Explicit(Rung-Kutta)/adaptive local time stepping      //!
!// 10-Steady/Unsteady: 	                    Unsteady                                    //!
!// 11-Gradient Calculation Scheme:	            non                                                    //!
!// 12-Turbulence Model:	                    non                                                    //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Inviscid_AdaptLocalTimStp
 Implicit None

!===============================
 Integer,Parameter::Dim=120000
!===============================
 Integer::I,J,K
 Integer::NS  !!  Number of Steps in rung-kutta method
 Real(8)::T1,T2 !This parameters show start and end time of code.
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
 Integer::UnStdy_cyc!!  Unsteady Cycle Steps
 Integer::M_max  !!  maximum Number of Internal steps
 Integer::PassNo
 Integer::NPassCell
 Integer::Cnt
 Integer::Tmp_NF !!Copy of Number of Faces of Mesh  
 Real(8)::GM   !Gama Constant (Specific Heat Ratio)
 Real(8)::Co   !!  Coefficients of residuals in Rung-kutta steps
 Real(8)::ALF  !Infinite Flow Angle to X Axis
 Real(8)::R0   !Density of Infinite Flow
 Real(8)::P0   !Pressure of Infinite Flow
 Real(8)::C0 !Sound Speed of Infinite Flow
 Real(8)::U0,V0   !Infinite Flow Velocity in X Direction
 Real(8)::E0 !Internal Energy of Infinite Flow
 Real(8)::CFLx !Currant Number for Explicit Methods
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::RM !Residual of Mass Equation
 Real(8)::U    !Velocity in x Direction
 Real(8)::V    !Velocity in Y Direction
 Real(8)::RKco !  Rung-kutta Coefficients
 Real(8)::TotTime  !Total Time for Unsteady simulation
 Real(8)::Sum_DT !!  The Sum of Time steps(Total time Passed For converging the Solution)
 Real(8)::UnStdyResMas !!Residual of Real Time Iteration
 Real(8)::DT_min !Minimum Time Step
 Real(8)::Temp!!  value of denominator in the friction of caculating Conservative Values at (N+1)th Time Step by dual time step method
 Real(8)::Interpolate2 !!  Local Vriable substituted Rung-kutta Coefficients
 Real(8)::InterpolateCoe !!Interpolate Coefficient in linear Interpolating
 Real(8)::ERMX !Maximum Error in Steady Approach
 Real(8)::RemainingTim !  Remaining Time Of the Total time considered for Solution 
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:Dim)::IPassCell
 Integer,Dimension(1:Dim)::NCELL_EDGE
 Integer,Dimension(1:Dim)::M !Number of characteristics
 Integer,Dimension(1:Dim)::Interpolate!!  Components Of linear Interpolation
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:Dim)::CELL_EDGE
 Integer,Dimension(1:Dim)::Tmp_IDS!!Copy of faces information 
 Real(8),Dimension(1:4,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:4,1:Dim)::WN !Conservative Values at (N)th Time Step
 Real(8),Dimension(1:4,1:Dim)::Wstor!!  Conservative Variables in Same Time Of Run
 Real(8),Dimension(1:4,1:Dim)::WC !Constant values Of Rung-kutta Method
 Real(8),Dimension(1:4,1:Dim)::Con !Convection Term of Mean flow Equations
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xc,Yc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::A !Area of each cell
 Real(8),Dimension(1:Dim)::NX,NY !Normal Vectors of each Face
 Real(8),Dimension(1:Dim)::DA   !Area of each Face
 Real(8),Dimension(1:Dim)::DT !Time step
 Real(8),Dimension(1:Dim)::P !Pressure
 Real(8),Dimension(1:Dim)::Pstor!!  Pressure Value of The Cells at the same time of the Run
 Real(8),Dimension(1:Dim)::Tmp_NX !Copy of Normal Vectors of each Face  
 Real(8),Dimension(1:Dim)::Tmp_NY!Copy of Normal Vectors of each Face  
 Real(8),Dimension(1:Dim)::Tmp_DA!Copy of Normal Vectors of each Face  
 Real(8),Dimension(1:5,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Integer,Dimension(1:10,1:2)::FaceNum!Copy of Index of Faces belongs to each region
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Integer,Dimension(1:4,1:Dim)::Corn !Corners point index of Each Cell 
 Real(8)::Rinf !Reynolds Number of infinite Flow 
 Real(8)::TT !Total Temprature
 Real(8)::K2,K4 !2nd and 4th derivative coefficient in Jameson artificial dissipation
 Real(8)::Time_Coe !Time Coefficient for dual time stepping
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient 
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
   
!part 8:
 Do j=1,NC
    
    IF(XC(j)<0.3)Then
      WNP1(1,J) = 1.0 
      WNP1(2,J) = 0.75
      WNP1(3,J) = 0.0
      WNP1(4,J) = 1.0/((GM-1))+ 0.5*(0.75*0.75 )
    Else
      WNP1(1,J) = 0.125 
      WNP1(2,J) = 0.0
      WNP1(3,J) = 0.0
      WNP1(4,J) = 0.1/((GM-1))
    Endif
    
 End do 
  
!Part 9: 
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do
 
!Part 10:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_RiemannForShockTube(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,WNP1,P,WB,XC) !Just For Shock Tube
 Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,WNP1,P,WB)
 Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,WNP1,P,WB)

!Part 11:
 Call CopyMesh(Dim,NF,IDS,NX,NY,DA,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
               Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum)
 
Call CPU_Time(t1)

!Part 12:
 Sum_DT=0
 UnStdy_cyc = 0

!Part 13:
 Do While( Sum_DT<TotTime )

   !Part 14:
	Call TimSTP_Inviscid(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,DT)

   !Part 15:
	DT_min = 1000.0
    Do J=1,NC
       IF( DT_min > DT(J) ) DT_min = DT(J)
	End Do
   	
   !Part 16:
    M_max=0
    Do J=1,NC
	   M(J)	=Int(Log10(DT(J)/DT_min)/Log10(2.0))
	   DT(J)=DT_min*2**(M(J))
	   IF( M_max < M(J) ) M_max =M(J)
    End Do

   !Part 17:
    RemainingTim = TotTime - Sum_DT

    IF( RemainingTim<DT_min*2**(M_max) )Then
	 DT_min = RemainingTim / (2**M_max)

     M_max=0
     Do J=1,NC
	    M(J) =Int(Log10(DT(J)/DT_min)/Log10(2.0))
	    DT(J)=DT_min*2**(M(J))
	    IF( M_max < M(J) ) M_max =M(J) 
     End Do
     
	 Print*,Sum_DT+DT_min*2**(M_max),'finished time'
    EndIF
  
   !Part 18:
    Do J=1,NC
	   WN(1,J)=WNP1(1,J)
	   WN(2,J)=WNP1(2,J)
	   WN(3,J)=WNP1(3,J)
	   WN(4,J)=WNP1(4,J)
    End Do
    
   !Part 19:
    Do PassNo=0, 2**(M_max)-1
	   NPassCell=0
	   Do J=1,NC
	      IF(MOD(PassNo,2**(M(J)))==0) then
	       NPassCell=NPassCell+1
	       IPassCell(NPassCell)=J
	      End IF !MOD 
	    End Do !J
		
	   !Part 20:
	    Call VectorizeMesh(Dim,NPassCell,IPassCell,NCELL_EDGE,CELL_EDGE,Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum,&
                       IDS,NX,NY,DA,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFIF1,NFIF2)
     
       !Part 21:
    	Do K=1,NPassCell
		   J=IPassCell(K)
		   WC(1,J) = WNP1(1,J)
           WC(2,J) = WNP1(2,J)
           WC(3,J) = WNP1(3,J)
           WC(4,J) = WNP1(4,J)	
		   Interpolate(J)=0    
        End Do

       !Part 22:
		Interpolate2=0

	   !Part 23:
	    Do NS=1,NRKS
           RKco=RKJ(NS)
		   cnt=0
	       Do J=1,NC
			  			 
			  K=IPassCell(cnt+1)
		      IF(J==K)Then
               Wstor(1,J) = WNP1(1,J)
               Wstor(2,J) = WNP1(2,J)
               Wstor(3,J) = WNP1(3,J)
               Wstor(4,J) = WNP1(4,J) 
               Pstor(J)   = P(J)  
			   cnt=cnt+1
		      Else
               InterpolateCoe=(Interpolate(J)+Interpolate2)*(DT_min/DT(J)) !RKco
		 	   Wstor(1,J)=InterpolateCoe*(WNP1(1,J)-WC(1,J))+WC(1,J)
			   Wstor(2,J)=InterpolateCoe*(WNP1(2,J)-WC(2,J))+WC(2,J)
			   Wstor(3,J)=InterpolateCoe*(WNP1(3,J)-WC(3,J))+WC(3,J)
			   Wstor(4,J)=InterpolateCoe*(WNP1(4,J)-WC(4,J))+WC(4,J)

               U    = Wstor(2,J)/Wstor(1,J)
               V    = Wstor(3,J)/Wstor(1,J)
               Pstor(J) = (GM-1)*(Wstor(4,J)-0.5*Wstor(1,J)*(U*U+V*V))
		      End IF 		  
              
           End Do
		   			   
		  !Part 24:
	       Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,Pstor,WB)
           Call BC_RiemannForShockTube(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,Wstor,Pstor,WB,XC) !Just For Shock Tube
           Call BC_InFlow(Dim,NFI1,NFI2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,Wstor,Pstor,ALF,Minf,WB)
           Call BC_VisOutFlow(Dim,NFO1,NFO2,IDS,GM,P0,Wstor,Pstor,WB)
           Call BC_Symmetry(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,Wstor,Pstor,WB)
          
           !part 25:
           Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,Wstor,WB,Pstor,Con)
          
		  !Part 26:
           Do K=1,NPassCell
		      J=IPassCell(K)	

	          Co =RKco*DT(J)/A(J)
		      WNP1(1,J) = WC(1,J) - Co* Con(1,J) 
		      WNP1(2,J) = WC(2,J) - Co* Con(2,J)
              WNP1(3,J) = WC(3,J) - Co* Con(3,J)
              WNP1(4,J) = WC(4,J) - Co* Con(4,J)
              
             !Part 27:
              U    = WNP1(2,J)/WNP1(1,J)
              V    = WNP1(3,J)/WNP1(1,J)
              P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
           End Do

		   !Part 28:
		   Interpolate2=RKco
		   
        End Do !Ns

		!Part 29:
		Do J=1,NC
		  Interpolate(J)=Interpolate(J)+1
		End Do
  	    	
	End Do !PassNo while

	!Part 30:
	Call UnVectorizeMesh(Dim,NF,IDS,NX,NY,DA,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,&
                     Tmp_IDS,Tmp_NX,Tmp_NY,Tmp_DA,Tmp_NF,FaceNum)

   !Part 31:
    UnStdy_cyc = UnStdy_cyc + 1
	Sum_DT = Sum_DT+ DT_min*2**(M_max)
    Print*,'UnStdy_cyc:',UnStdy_cyc,'  Time:',Sum_DT

   !Part 32:
	Do J=1,NC
        DT(J) = DT_min*2**(M_max)
    End Do
    
    !Part 33:
    Call ResMass(Dim,NC,WNP1,WN,DT,Ncyc,Rm)
		
   !Part 34:
    IF( Mod(UnStdy_cyc,NWrite)==0 )Then
     Print*,'Writing Results... '
     Call Write_ShockTubeResults(Dim,NC,Xc,WNP1,P) !Just For Shock Tube
	End If
   
 End Do !Unsteady

!Part 35:
 Call Write_ShockTubeResults(Dim,NC,Xc,WNP1,P) !Just For Shock Tube
 
 Call CPU_Time(t2)
 Print*,'CPU Time: ' , (t2-t1)
!**********************************************************************************************************
 End 
!##########################################################################################################









