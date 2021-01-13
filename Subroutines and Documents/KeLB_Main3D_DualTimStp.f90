!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Feb., 20, 2016                                                                   //!
!// Developed by: *//*-+/                       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeLB_Main3D_DualTimStp(Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                      NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Wb,WNP1,Mu,P,GM,DUY,WTNM1,WTN,DT_Real,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,X,Y,Z,NX,NY,NZ,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                NFF1,NFF2,NP,IDS,FaceType,INW,XC,YC,ZC,DW,DT,Vol,MR,NRKS,RKJ,Wb,WNP1,Mu,P,GM,DUY,WTNM1,WTN,DT_Real
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,I,II,J,NS,ME
 Real(8)::RKco,AA,Temp,K,Epsilon,Rho,Co,tauwall,ustar,yplus,Yn,fmu,DT_Real
 Real(8)::Ceps1,Ceps2,SigK,Sige,CMu,Lam1,Lam2,Lam3
!-------------------------------------------------------------------------------------------
 Integer::NC
 Integer::NP
 Integer::NF !Number of Faces Constructing Mesh
 Integer::NRKS
 Integer::NFW1,NFW2
 Integer::NF1,NF2
 Integer::NFI1,NFI2
 Integer::NFO1,NFO2
 Integer::NFS1,NFS2
 Integer::NFF1,NFF2
 Real(8)::MR
 Real(8)::GM
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTNM1,WTN,Rest
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:Dim)::Vol
 Real(8),Dimension(1:Dim)::Mu
 Real(8),Dimension(1:Dim)::Mut
 Real(8),Dimension(1:Dim)::DT
 Real(8),Dimension(1:2,1:Dim)::WTB
 Real(8),Dimension(1:2,1:Dim)::WTC
 Real(8),Dimension(1:2,1:Dim)::Cont
 Real(8),Dimension(1:2,1:Dim)::Dift
 Real(8),Dimension(1:2,1:Dim)::St
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:Dim)::NX,NY,NZ
 Real(8),Dimension(1:Dim)::XC,YC,ZC
 Real(8),Dimension(1:Dim)::DW
 Real(8),Dimension(1:Dim)::DUY
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
 Real(8),Dimension(1:Dim)::DKX_F,DKY_F,DKZ_F
 Real(8),Dimension(1:Dim)::DEPSX_F,DEPSY_F,DEPSZ_F
 Real(8),Dimension(1:5)::RKJ
!*********************************************************************************************
!Part 1:
 Ceps1    =1.44D0
 Ceps2    =1.92D0
 SigK     =1.0D0
 Sige     =1.3D0
 CMu      =0.115D0
 Lam1     =1.0D0
 Lam2     =0.4D0
 Lam3     =0.2D0

!Part 2:
 Do I=1,NC
    WTC(1,I) =WTNP1(1,I)
    WTC(2,I) =WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
      
   !Part 4:    
	RKco=RKJ(NS)
       
   !Part 5:
    Call KeLam_BC3D(Dim,NX,NY,NZ,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB)

   !Part 6:  
    Call Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
         
   !Part 7:  
    Call KFi_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1,WTNP1,WB,WTB,DKX_F,DKY_F,DKZ_F,DEPSX_F,DEPSY_F,DEPSZ_F)  

   !Part 8:
	Call KFi_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)

   !Part 9:
	Call KFi_Dif3d(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DEPSX_F,DEPSY_F,DEPSZ_F,Sigk,Sige,MR,Mu,Mut,Dift)
         
   !Part 10:
    Call KeLB_Source3D(Dim,NC,IDS,DW,Vol,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,&
                     GM,P,WTNP1,Mu,Mut,DUY,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,St)

   !Part 11:
    Do J=1,NC
        
        AA = Vol(J)
        Rest(1,J) = -(Cont(1,J)+Dift(1,J)+St(1,J))/AA
        Rest(2,J) = -(Cont(2,J)+Dift(2,J)+St(2,J))/AA
        
        Temp = 2*DT_Real+3*RKco*DT(J)
        
        WTNP1(1,J) = ( 2*WTC(1,J)*DT_Real + RKco*DT(J)*(4*WTN(1,J) - WTNM1(1,J) + 2*Rest(1,J)*DT_Real) ) / Temp
        WTNP1(2,J) = ( 2*WTC(2,J)*DT_Real + RKco*DT(J)*(4*WTN(2,J) - WTNM1(2,J) + 2*Rest(2,J)*DT_Real) ) / Temp
        
      !Part 12: 
       if(WTNP1(1,J)<0.0 )   WTNP1(1,J)=WTC(1,J)    !        WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 )   WTNP1(2,J)=WTC(2,J)    !        WTNP1(2,J)=WTN(2,J)
	   
      !Part 13:
       Rho     = WNP1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho

      !Part 14:
       II = INW(J)     ! II: Wall Face
       Yn = DW(J) 
       ME = IDS(1,II)

       Tauwall = Mu(ME)*DUY(II) !!!!!!!!!!!!!!!!!!!!
       Ustar   = Dsqrt(abs(MR*Tauwall/WNP1(1,ME)))
       Yplus   = (1.0/MR)*Rho*Ustar*Yn/Mu(J)

	  
	   FMu = 0.04 + (1.0 - 0.04)*((1.0 - EXP(-(Yplus - 8.0)/26.0))**2)
	   
      !Part 15:  
       Mut(J) =(1.0/MR)*Cmu*FMU*Rho*K*K / (Epsilon+1.e-20)
	  
      End Do !J  
 End Do !NS
!*********************************************************************************************
 End
!###########################################################################################
