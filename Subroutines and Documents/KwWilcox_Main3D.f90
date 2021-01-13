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
!// Date: May., 15, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KwWilcox_Main3D(Dim,NC,NP,INW,IDS,MR,NRKS,RKJ,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,DA,FaceType,&
                            NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,Z,NX,NY,NZ,XC,YC,ZC,DW,Vol,Mu,DT,Mut)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NC,NP,INW,IDS,MR,NRKS,RKJ,P,WNP1,WB,GM,NF,NF1,NF2,DA,FaceType,NFW1,NFW2,&
                 NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,Z,NX,NY,NZ,XC,YC,ZC,DW,Vol,Mu,DT
 Intent (Out  )::WTNP1,Mut
 
 Integer::I,J,ME,NS
 Real(8)::RKco,CO,Omega,K,Rho,ALFA,BETA,BETA_S,Sigk,Sigw
!-------------------------------------------------------------------------------------------
Integer::Dim
Integer::NC
Integer::NF !Number of Faces Constructing Mesh
Integer::NF1,NF2
Integer::NRKS
Integer::NP
Integer::NFW1,NFW2
Integer::NFF1,NFF2
Integer::NFI1,NFI2
Integer::NFO1,NFO2
Integer::NFS1,NFS2
Real(8)::MR
Real(8)::GM
Integer,Dimension(1:6,1:Dim)::IDS
Integer,Dimension(1:Dim)::FaceType
Integer,Dimension(1:Dim)::INW
Real(8),Dimension(1:2,1:Dim)::WTN
Real(8),Dimension(1:2,1:Dim)::WTNP1
Real(8),Dimension(1:2,1:Dim)::WTB
Real(8),Dimension(1:2,1:Dim)::Cont
Real(8),Dimension(1:2,1:Dim)::Dift
Real(8),Dimension(1:2,1:Dim)::St
Real(8),Dimension(1:5,1:Dim)::WNP1
Real(8),Dimension(1:6,1:Dim)::WB
Real(8),Dimension(1:Dim)::Mut
Real(8),Dimension(1:Dim)::P
Real(8),Dimension(1:Dim)::DA
Real(8),Dimension(1:Dim)::X,Y,Z
Real(8),Dimension(1:Dim)::NX,NY,NZ
Real(8),Dimension(1:Dim)::XC,YC,ZC
Real(8),Dimension(1:Dim)::DW
Real(8),Dimension(1:Dim)::Vol
Real(8),Dimension(1:Dim)::Mu
Real(8),Dimension(1:Dim)::DT
Real(8),Dimension(1:5)::RKJ
Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C
Real(8),Dimension(1:Dim)::DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,DOmegZ_F
!***************************************** Turb_Main ********************************************
!Part 1:
 ALFA   = 5.0/9.0
 BETA   = 3.0/40.0
 BETA_S = 0.09
 Sigk   = 0.5
 Sigw   = 0.5

!Part 2:
 Do I=1,NC
    WTN(1,I) = WTNP1(1,I)
    WTN(2,I) = WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
      
   !Part 4:    
	RKco = RKJ(NS)
       
   !Part 5:
    Call Kw_BC3D(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,NX,NY,NZ,DW,Mu,WB,WNP1,WTNP1,WTB)

   !Part 6:  
    Call Velocity_CellGrad3D(Dim,NC,NF,NF1,NF2,IDS,Vol,NX,NY,NZ,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C)
            
   !Part 7:  
    Call KFi_FaceGrad3D(Dim,NFW1,NFW2,NF,NF1,NF2,NFS1,NFS2,NP,NC,IDS,FaceType,X,Y,Z,XC,YC,ZC,WNP1,WTNP1,WB,WTB,DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,DOmegZ_F)

   !Part 8:
    Call KFi_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)
  
   !Part 9:
    Call KFi_Dif3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,DOmegZ_F,Sigk,Sigw,MR,Mu,Mut,Dift)

   !Part 10:
    Call KwWilcox_Source3D(Dim,NC,MR,ALFA,BETA,BETA_S,Mut,Vol,WTNP1,WNP1,WB,DUX_C,DUY_C,DUZ_C,DVX_C,DVY_C,DVZ_C,DWX_C,DWY_C,DWZ_C,St)
      
   !Part 11:
    Do J=1,NC

	   Co = RKco*DT(J)/Vol(J)
        
       WTNP1(1,J) = WTN(1,J)  - Co*( ConT(1,J) + DifT(1,J) - St(1,J) ) 
       WTNP1(2,J) = WTN(2,J)  - Co*( ConT(2,J) + DifT(2,J) - St(2,J) )
	    
      !Part 12: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)

      !Part 13:
       Rho       = WNP1(1,J)
       K         = WTNP1(1,J)/Rho
       Omega     = WTNP1(2,J)/Rho
         
      !Part 14:  
       Mut(J) =(1.0/MR)*Rho*K / (Omega+1.e-20)

   End Do !J

 End Do !NS
  
!*********************************************************************************************
 End
!###########################################################################################