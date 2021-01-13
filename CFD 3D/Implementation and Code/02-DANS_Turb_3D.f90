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
!// gradients. This code uses wide range wide range of turbulence models. ranging for calculating      //!
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
!// 9-Transient Term Discretization Scheme:     Explicit(Rung-Kutta)                                   //!
!// 10-Steady/Unsteady: 	                    Steady and Unsteady                                    //!
!// 11-Gradient Calculation Scheme:	            Green-Gauss                                            //!
!//                                             Least Square                                           //!
!// 12-Turbulence Model:                        KWSST                                                  //!
!//                                             KWSST                                                  //!
!//                                             KeLB                                                   //!
!//                                             KwWilcox                                               //!
!//                                             LES_WALE                                               //!
!//                                             LES_DSmag                                              //!
!//                                             KWSST_Transition                                       //!
!//                                             SA                                                     //!
!// 13-Moving Mesh Method:	                    non                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program DANS_Turb_3D
 use MainVariables
 use MeshBC3DHeader
 use Read_3DMeshHeader
 use FaceOfCellHeader
 Implicit None

 Integer::I,J,Dim
 nLocalArrayFace=10
 nLocalArrayCell=10
!***************************************** Main ********************************************

 Call MeshVarConstructir("mesh.gid")
 
 SumCF=0.0;SumCP=0.0;sum_umean=0.0;sum_vmean=0.0;sum_uvmean=0.0;sum_u2mean=0.0;sum_v2mean=0.0
 
 
!Part 4:
 Call Read_3DMesh(NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z,"mesh.gid")
      
!Part 5:
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"Imesh.plt")
 print*,'read mesh finished'
 
!Part 6:
 Call Read_Setting(Minf,Rinf,ALF,Tt,ERmx,CFLx,NRKS,NWrite,Init,K2,K4,RKJ,TotTime,Time_Coe)

!Part 7:
 Call MeshBC3D(IntegerLocalArrayFace,NR,NFR,BC,IDS,FaceType,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)

!Part 8:
 Call FaceOfCell(NF,NC,IDS,NFace_Cell,IFace_Cell)
 
!Part 9:
 DO cell=1,NC
    CALL PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 END DO
 print*,'find points of cell finished'
 
!Part 10:
 Call GeoCal3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Vol,DA,Nx,Ny,Nz,XC,YC,ZC)

!Part 11:
 Call WallDistance3D(Dim,NC,NFW1,NFW2,FaceType,IDS,X,Y,Z,Xc,Yc,Zc,INW,DW)
 print*,'pre proc finished'
 
 
!Part 12:
 Call InitMeanFlow3D(Dim,Init,NC,ALF,Minf,Rinf,Tt,MR,GM,PrL,PrT,R0,P0,C0,U0,V0,W0,T0,Mu0,B0,WNP1)

!Part 13:
 Do J=1,NC 
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    W = WNP1(4,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))
    
   !Part 14:
    Temp = GM*P(J)/WNP1(1,J) 

   !Part 15:
    Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)     

 End Do

!Part 16:
 Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
 Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
 Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)

!Part 17:
 Ncyc = 0
 Rm   = 10.0
 Naverage = 0
 Time = 0.0
 
!Part 18:
 !!!Call Kw_Init(Dim,NC,MR,WTNP1,Mut)
 !!!Call Ke_Init(Dim,NC,MR,WTNP1,Mut) 
 Call KwSST_Trans_Init3D(Dim,NC,IDS,Minf,Rinf,Wtnp1,Mut,TUinf,Mut0,Gamainf,Kinf,Omegainf,WNP1) 
 !!!Call SA_Init3D(Dim,NC,Mu0,R0,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mut0,Nuhat0,WTNP1,Mut)
 !!!TurbQ(:) =  0.0
 !!!Mut(:)=0.0
 
!Part 19:
 Do While(Rm > ERmx)

   !Part 20:
    Ncyc=Ncyc+1
     
   !Part 21:
    Do J=1,NC
       WC(1,J) = WNP1(1,J)
       WC(2,J) = WNP1(2,J)
       WC(3,J) = WNP1(3,J)
       WC(4,J) = WNP1(4,J)
       WC(5,J) = WNP1(5,J)   
    End Do

   !Part 22:
	Call TimSTP_Turb3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,VOL,CFLx,GM,P,WNP1,WB,Mu,Mut,PrL,PrT,MR,DT)

   !Part 23:
    !!!DT_min = 1000.0
    !!!Do J=1,NC
    !!!   IF( DT_min > DT(J) ) DT_min = DT(J)
    !!!End Do
    !!!DT(:) = DT_min
    !!!Time = Time + DT_min

   !Part 24:
    Do NS=1,NRKS
   
      !Part 25:
	   RKco=RKJ(NS)

      !Part 26:
       !!!Call FirstOrd_Gradient3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,WNP1,WB,P,GWNP1)
       !!!Call Limiter3D(Dim,NC,NF1,NF2,IDS,GM,XC,YC,ZC,FaceType,WNP1,P,WB,GWNP1,Limit,X,Y,Z,NF)
       !!!Call ConMeanFlow_AUSM_Plus_HO3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,FaceType,DA,GM,WNP1,WB,P,GWNP1,Xc,Yc,Zc,Limit,Con,X,Y,Z)
       
      !Part 27:
	   Call ConMeanFlow_AUSM_PlusUP3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,Minf,WNP1,WB,P,Con)
       !!!Call ConMeanFlow_ScalarDiss_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,Con)
       
      !Part 28:
       Call VelTemp_GradFace3D(Dim,NC,NF,NP,NF1,NF2,NFW1,NFW2,NFS1,NFS2,IDS,X,Y,Z,Xc,Yc,Zc,FaceType,&
                              WNP1,WB,GM,P,NX,NY,NZ,&
                              DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ,DTX,DTY,DTZ)
    
      !Part 29:
       Call DifMeanFlowTurbNoWallFu3D(Dim,NC,NF1,NF2,NFW1,NFW2,NF,IDS,GM,PrL,PrT,MR,Mu,Mut,NX,NY,NZ,DTX,DTY,DTZ,DUX,DUY,DUZ,&
                                      DVX,DVY,DVZ,DWX,DWY,DWZ,WNP1,WB,TurbQ,Dif)
            
      !Part 30:
       Do J=1,NC

		  Co = RKco*DT(J)/Vol(J)
         
          WNP1(1,J) = WC(1,J) - Co* ( Con(1,J)            )
		  WNP1(2,J) = WC(2,J) - Co* ( Con(2,J) + Dif(2,J) )
          WNP1(3,J) = WC(3,J) - Co* ( Con(3,J) + Dif(3,J) )
          WNP1(4,J) = WC(4,J) - Co* ( Con(4,J) + Dif(4,J) )
          WNP1(5,J) = WC(5,J) - Co* ( Con(5,J) + Dif(5,J) )
          
         !Part 31:
	      U = WNP1(2,J)/WNP1(1,J)
          V = WNP1(3,J)/WNP1(1,J)
          W = WNP1(4,J)/WNP1(1,J)
          P(J) = (GM-1)*(WNP1(5,J)-0.5*WNP1(1,J)*(U*U+V*V+W*W))
          
         !Part 32:
          Temp = GM*P(J)/WNP1(1,J) 
          
         !Part 33:
          Mu(j) = (Temp**1.5)*(1.0+B0)/(Temp+B0)     

       End Do

      !Part 34: 
       Call BC_Wall3D(Dim,NFW1,NFW2,IDS,GM,P,WB)
       Call BC_Riemann3D(Dim,NFF1,NFF2,GM,U0,V0,W0,P0,R0,C0,IDS,WNP1,Nx,Ny,Nz,DA,P,WB)
       Call BC_Symmetry3D(Dim,NFS1,NFS2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
       Call BC_InFlow3D(Dim,NFI1,NFI2,NX,NY,NZ,DA,IDS,GM,U0,V0,W0,P0,R0,WNP1,P,ALF,Minf,WB)
       Call BC_VisOutFlow3D(Dim,NFO1,NFO2,NX,NY,NZ,DA,IDS,GM,P0,WNP1,P,WB)

    End Do !Ns	

   !Part 35:
   !=====================================================================================================================
    !!!Call KWSST_Main3D(Dim,Ncyc,INW,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
    !!!                  NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call KeLB_Main3D(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
    !!!               NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Wb,WNP1,Mu,P,GM,DUY,WTNP1,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call KwWilcox_Main3D(Dim,NC,NP,INW,IDS,MR,NRKS,RKJ,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,DA,FaceType,&
    !!!                     NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,Z,NX,NY,NZ,XC,YC,ZC,DW,Vol,Mu,DT,Mut)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call LES_WALE_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,DW,Mut)
    !!!TurbQ(:) =  0.0

    !!!Call LES_DSmag_Main(Dim,NC,NF1,NF2,NF,MR,Vol,IDS,WNP1,WB,NX,NY,NZ,Mut,TurbQ)
    !!!
    !!!Call KWSST_Transition_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,&
    !!!                             NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,&
    !!!                             DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut,Kinf,Omegainf,Gamainf)
    !!!TurbQ(:) =  WTNP1(1,:)
    
    !!!Call SA_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
    !!!               NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,WB,WNP1,&
    !!!               Mu0,Mu,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Nuhat0,WTNP1,Mut)
    !!!TurbQ(:) =  0.0
    
    Call KWSST_Transition_Main3D(Dim,Ncyc,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,&
                                 NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,FaceType,XC,YC,ZC,DW,&
                                 DT,Vol,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut,Kinf,Omegainf,Gamainf)
   !=====================================================================================================================

    !Part 36:
 	 Call ResMass3D(Dim,NC,WNP1,WC,DT,Ncyc,Rm)
     Print*,'Ncyc:',Ncyc,'Mass eq. Residual:',Rm,maxval(Mut )

    !Part 37:
     IF( Mod(Ncyc,NWrite)==0 )Then
      Print*,'Writing Results... ',Ncyc,Rm

      Call Write_CP3D(Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P)
      Call Write_CF3D(Dim,Minf,Rinf,NFW1,NFW2,X,Y,Z,IDS,DUY,Mu)
      
      Call Write_FlatPlateAnalyticData3D(Dim,NFW1,NFW2,NC,IDS,X,Xc,Yc,Minf,Rinf,WNP1,Mu,DUY)
      
      Call Write_VelocityContour3D(Dim,NC,NP,X,Y,Z,Corn,WNP1)
      
      Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,P,"Pressr.plt")
      Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,Mut,"Muturb.plt")
      
      Call Write_ConservativeVariables3D(Dim,NC,WNP1)
     End If
     
 End Do !Do While
 
!*********************************************************************************************
 End 
!###########################################################################################


