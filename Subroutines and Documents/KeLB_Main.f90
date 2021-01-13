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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: H. Nazari, Mechanical Eng., Amirkabir University of Technology           //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeLB_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                      NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,P,GM,DUY,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,&
                NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,P,GM
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,I,II,J,NC,NS,NP,NRKS,ME,NF,NFW1,NFW2,NF1,NF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::RKco,K,Epsilon,Rho,MR,Co,tauwall,ustar,yplus,Yn,fmu
 Real(8)::Ceps1,Ceps2,SigK,Sige,CMu,Lam1,Lam2,Lam3,GM
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,WTN,Cont,Dift,St
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,XC,YC,DW,A,Mu,Mut,DT,DUY,DUX_C,DUY_C,DVX_C,DVY_C,&
                           DKX_F,DKY_F,DEPSX_F,DEPSY_F,P
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
    WTN(1,I) =WTNP1(1,I)
    WTN(2,I) =WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
      
   !Part 4:    
	RKco=1.0/(NRKS-NS+1)
       
   !Part 5:
    Call KeLam_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB)

   !Part 6:  
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
         
   !Part 7:  
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX_F,DKY_F,DEPSX_F,DEPSY_F) 
	 
    Do I=NFS1+1,NFS2
       DKY_F(I)=0.0
	   DEPSY_F(I)=0.0
    End do

   !Part 8:
	Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,WTNP1,Wb,WTB,Cont)

   !Part 9:
	Call KFi_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,MR,Mu,Mut,Dift)
         
   !Part 10:
    Call KeLB_Source(Dim,NC,IDS,DW,A,INW,MR,Ceps1,Ceps2,WNP1,Lam1,Lam2,Lam3,&
                     GM,P,WTNP1,Mu,Mut,DUY,DUX_C,DUY_C,DVX_C,DVY_C,St)

   !Part 11:
    Do J=1,NC
        
	   Co = RKco*DT(J)/A(J)
         
       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) + St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) + St(2,J) ) 
    
      !Part 12: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)
	   
      !Part 13:
       Rho     = WNP1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho

      !Part 14:
       II = INW(J)     ! II: Wall Face
       Yn = DW(J) 
       ME = IDS(1,II)

       Tauwall = Mu(ME)*DUY(II)
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
