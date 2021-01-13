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
!// Cheif Developer: *+-+-*, Amirkabir University of Technology                                //!
!// *-+*                                                                                               //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.FFFF.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!//                                                                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program Airfoil_Optimization_Plasma_Main
 Implicit None
!******************************************************************************************* 
 Integer,Parameter::Dim = 60000
 Integer,Parameter::SwarmSize = 5		! Number of Swarms considered 
 Integer,Parameter::Iterations =300	    ! Particle swarm optimization Maximum iterations number
 Real(8),Parameter::ConvergERROR = 1e-2	! Particle swarm optimization Convergance condition
 Integer,Parameter::NVariable=2
 
 Integer::I,J,S,K,P1
 Real(8)::T0 !Temperature of Infinite Flow
 Real(8)::E0
 Real(8)::B0
 Real(8)::Rinf
 Real(8)::Tt
 Real(8)::MR
 Real(8)::PrL
 Real(8)::PrT
 Real(8)::RKco
 Real(8)::Mu0
 Real(8)::Mut0
 Real(8)::K2,K4
 Real(8)::time
 Real(8)::Temp
 Real(8)::TotTime
 Real(8)::ERmxUnStdy
 Real(8)::Time_Coe    
 Real(8),Dimension(1:4,1:Dim)::Dif ! Diffusion Term of Mean flow Equations
 Real(8),Dimension(1:Dim)::DW  !Distance to Nearest Wall
 Real(8),Dimension(1:Dim)::Mu
 Real(8),Dimension(1:Dim)::Mut
 Real(8),Dimension(1:Dim)::DUX,DUY,DVX,DVY,DTX,DTY
 Real(8),Dimension(1:Dim)::Taukk
 Real(8),Dimension(1:Dim)::DT
 Integer,Dimension(1:Dim)::INW  !Index of Nearest Wall
 Real(8),Dimension(1:2,1:Dim)::WTNP1 !Turbulence Variables at new time step
 Real(8),Dimension(1:5)::RKJ !Runge Kutta Jameson Coefficient 
 Integer::NP
 Integer::NC
 Integer::NF
 Integer::NR
 Integer::NBP
 Integer::UpRegion
 Integer::LwRegion
 Integer::DimU,DimL
 Integer::BSOrder_Up
 Integer::BSOrder_Lw
 Integer::NPtCurv_Up
 Integer::NPtCurv_Lw
 Integer::Check
 Integer::NF1,NF2
 Integer::NFW1,NFW2
 Integer::NFF1,NFF2
 Integer::NFI1,NFI2
 Integer::NFS1,NFS2
 Integer::NFO1,NFO2
 Integer::NFIF1,NFIF2
 Integer::NRKS
 Integer::NWrite
 Integer::Init
 Integer::Ncyc
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points													  
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Integer,Dimension(1:4,1:Dim)::Corn !Corners point index of Each Cell 
 Integer,Allocatable,Dimension(:)::IBP_Lw,IBP_Up
 Real(8)::Zita_TE
 Real(8)::TE_Thick
 Real(8)::Shap
 Real(8)::Xo,Yo
 Real(8)::GM
 Real(8)::ALF
 Real(8)::R0
 Real(8)::P0
 Real(8)::C0
 Real(8)::U0
 Real(8)::V0
 Real(8)::ERmx
 Real(8)::CFLx
 Real(8)::Minf
 Real(8)::U,V
 Real(8)::CL
 Real(8)::CD
 Real(8)::ObjGbest
 Real(8)::ObjectiveFun
 Real(8),Dimension(1:Dim)::Xc,Yc !Coordinate of Center of Element
 Real(8),Dimension(1:Dim)::A
 Real(8),Dimension(1:Dim)::NX,NY
 Real(8),Dimension(1:Dim)::DA
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:Dim)::DelX,DelY
 Real(8),Dimension(1:4,1:Dim)::WNP1 !Conservative Values at (N+1)th Time Step
 Real(8),Dimension(1:5,1:Dim)::WB !Conservative Values(1:4,:) and Pressure(5:5,:) at Boundary Faces
 Real(8),Allocatable,Dimension(:)::Pol_Coeff_Up,Pol_Coeff_Lw
 Real(8),Allocatable,Dimension(:)::X_Up,Y_Up
 Real(8),Allocatable,Dimension(:)::X_Lw,Y_Lw
 Real(8),Allocatable,Dimension(:)::S_Up,S_Lw
 Real(8),Allocatable,Dimension(:)::Ynew_Lw,Ynew_Up
 Real(8),Allocatable,Dimension(:)::ShapeFunc_Coe
 Real(8),Allocatable,Dimension(:)::gBest
 Real(8),Allocatable,Dimension(:)::GeoCoeff
 Real(8),Allocatable,Dimension(:)::Pmax,Pmin
 Real(8),Allocatable,Dimension(:)::ObjPbest
 Real(8),Allocatable,Dimension(:)::ObjFuncS
 Real(8),Allocatable,Dimension(:,:)::BMI_Up,BMI_Lw
 Real(8),Allocatable,Dimension(:,:)::Pos
 Real(8),Allocatable,Dimension(:,:)::Velocity
 Real(8),Allocatable,Dimension(:,:)::pBest
 Real(8),Allocatable,Dimension(:,:)::CheckMatrix  
 Logical:: Converged = .False.
 Integer,Dimension(1:Dim)::Delta
 Real(8),Allocatable,Dimension(:)::PlasmaParameter
 Real(8)::freq
 Real(8)::PlasmaHeight
 Real(8)::Plasmawidth
 Real(8)::appliedVoltage
 Real(8)::PlasmaGap
 Real(8)::Roh_c
 Real(8)::e_c
 Real(8)::Eb
 Real(8)::DischargeTime
 Real(8),Dimension(1:Dim)::F_DBD_x,F_DBD_y
!*******************************************************************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)

!Part 3:
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 4:
 Call GeoCal2D(Dim,NF1,NF2,NF,NC,IDS,X,Y,Xc,Yc,NX,NY,DA,A)

!Part 5:
 Call InitMeanFlow(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,T0,Mu0,B0,WNP1)

!Part 6:
 Do J=1,NC
    U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
    Temp = GM*P(J)/WNP1(1,J) 
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)         
 End Do
 
!Part 7:
 Call Kw_Init(Dim,NC,MR,WTNP1,Mut)
 
!Part 8: 
 Call  WallDist(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,INW,DW)
 
!Part 9:
 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
!Part 10:
 Call PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,IDS,Corn)
 
!Part 11:
 Call PlasmaShayyParameters(freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,Roh_c,e_c,Eb,DischargeTime)

!Part 12:
Allocate( Pos(1:SwarmSize,1:NVariable),Velocity(1:SwarmSize,1:NVariable),&
		   pBest(1:SwarmSize,1:NVariable),gBest(1:NVariable),&
		   Pmax(1:NVariable),Pmin(1:NVariable),ObjpBest(1:SwarmSize),&
		   CheckMatrix(1:SwarmSize,1:NVariable),ObjFuncS(1:SwarmSize),PlasmaParameter(1:NVariable) )

 Open(16,File='objectiveFunctions.txt')
 Open(55,File='Converge Check.txt')

!Part 13:  
 Call PSO_Initial_plasma(Swarmsize,NVariable,Pmax,Pmin,Pos,Velocity,pBest)

!Part 14:
 Do I=1,Swarmsize
	ObjPbest(I) = 10000000 !minmizing
   !ObjPbest(I) = 0.0000001 !Maximizing 
 End Do

 ObjGbest = 10000000 !minimizing
!ObjGbest = 0.00000001 !Maximizing 

!Part 15:
 K = 0				
 Do While( K < Iterations .and. .Not.Converged  )
    K = K+1
    Print*,'Optimization Iterations:',K 

   !Part 16:
	Do S=1,SwarmSize
        
       Do I=1,NVariable
	      PlasmaParameter(I) = Pos(S,I)
       End Do

      !Part 17:
       Call PlasmaAfectedCell(Dim,PlasmaHeight,Plasmawidth,NC,Xc,Yc,ALF,BSOrder_Up,BSOrder_Lw,ShapeFunc_Coe,PlasmaParameter(1),Zita_te,Delta)
      
      !Part 18:
       Call PlasmaShayy(Dim,NC,freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,&
                        Roh_c,e_c,Eb,DischargeTime,Xc,Yc,Delta, F_DBD_x,F_DBD_y)
       
      !Part 19:
       Call Solver_Turb_withSource(&
             Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,&
             NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,Xc,Yc,NX,NY,DA,A,Dw,INW,&
             ERmx,CFLx,NRKS,NWrite,RKJ,&
             Minf,Rinf,MR,ALF,GM,R0,P0,T0,B0,C0,U0,V0,Tt,PrL,PrT,Mu0,Mut0, &
             WB ,WNP1,WTNP1,P,Mu,Mut,DUY,F_DBD_X,F_DBD_Y)
         
	  !Part 20: 
       Call PressLiftDragCo_Viscous(Dim,NFW1,NFW2,GM,Rinf,Minf,ALF,IDS,WNP1,WB,NX,NY,DA,Mu,DUY,CL,CD)
 
       ObjectiveFun = abs(Cd/Cl)

       Print*,' Swarm :' , S ,'ObjectiveFun' , ObjectiveFun 	 

	   ObjFuncS(S) =  ObjectiveFun

	  !Part 21:    
	   If( ObjectiveFun < ObjPbest(S) )Then
           
	    Do I=1,NVariable
           pBest(S,I) = Pos(S,I)
        End do
 
   		ObjPbest(S) = ObjectiveFun
       End If

      !Part 22: !For Maximizing use ">" !For minimizing use "<" 
      If(ObjPbest(S) <  ObjGbest )Then

	   Do I=1,NVariable
          gBest(I) = pBest(S,I)
       End do
        
       ObjGbest = ObjPbest(S)
       
      End If

    End Do

   !Part 23:
    Call PSO_Update_Plasma(Iterations,Swarmsize,NVariable,Pmax,Pmin,pBest,gBest,ObjFuncS,ObjGbest,Pos,Velocity)

   !Part 24:
 	Check = 0
    Do S=1,SwarmSize
       If (ABS(ObjFuncS(S) - ObjGbest)<= ConvergERROR) Check = Check + 1
    End Do
    If(Check == SwarmSize ) converged = .True.

   !Part 25:
    Write(55,'(I4,F10.6,F10.6,F10.6,F10.6,F10.6)') K , ObjFuncS(1) , ObjFuncS(2), ObjFuncS(3), ObjFuncS(4), ObjFuncS(5)

End Do !While
    
  
!Part 26:
 Do I=1,NVariable
	PlasmaParameter(I) = gBest(I)
 End Do
 
 Call PlasmaAfectedCell(Dim,PlasmaHeight,Plasmawidth,NC,Xc,Yc,ALF,BSOrder_Up,BSOrder_Lw,ShapeFunc_Coe,PlasmaParameter(1),Zita_te,Delta)
    
 Call PlasmaShayy(Dim,NC,freq,PlasmaHeight,Plasmawidth,appliedVoltage,PlasmaGap,&
                  Roh_c,e_c,Eb,DischargeTime,Xc,Yc,Delta, F_DBD_x,F_DBD_y)
       
 Call Solver_Turb_withSource(&
 Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,&
 NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2,Xc,Yc,NX,NY,DA,A,Dw,INW,&
 ERmx,CFLx,NRKS,NWrite,RKJ,&
 Minf,Rinf,MR,ALF,GM,R0,P0,T0,B0,C0,U0,V0,Tt,PrL,PrT,Mu0,Mut0, &
 WB ,WNP1,WTNP1,P,Mu,Mut,DUY,F_DBD_x,F_DBD_Y)

 Write(16,'(a)') 'final CL'							
 Write(16,'(F10.6)') CL 
 Write(16,'(a)') 'final Cd'							
 Write(16,'(F10.6)') Cd
 Write(16,'(a)') 'Final ObjGbest' 
 Write(16,'(F10.6)') ObjGbest
 
 Call Write_CF(Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu)
 Call Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
 Call Write_Contours(Dim,NC,NP,Corn,X,Y,GM,WNP1,P)

!*******************************************************************************************
End
!###########################################################################################
