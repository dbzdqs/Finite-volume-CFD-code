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
!// Supervisor:                                                                            //!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Date: Feb., 10, 2018                                                                   //!
!// Supervisor:                                                                            //!
!// Developed by: M. H. Saadat, Sharif University of Technology                            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KWSST_Main(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,&
                     NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,RKJ,Mu0,Wb,WNP1,Mu,WTNP1,Mut)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,MR,NRKS,RKJ,Mu0,WB,WNP1,Mu
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,I,J,Ii,Jj,NC,NF,NFW1,NFW2,NF1,NF2,NS,Ncyc,NP,NRKS,ME,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::Mu0,RKco,K,Omega,MR,Co,Vor,Sigk1,Sigk2,Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,Rho
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,WTN,Cont,Dift,St
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C,DKX_F,DKY_F,&
                           DOmegX_F,DOmegY_F,F11,F22,Sigk,Sigw,Beta,Gama,X,Y,NX,NY,A,Mu,Mut,XC,YC,DW,DT
 Real(8),Dimension(1:5)::RKJ
!*********************************************************************************************
!Part 1:
 Bstar=0.09
 Sigk1=0.85   ;   Sigw1=0.5     ;   Beta1=0.075    ;   Gama1=5.0/9.0
 Sigk2=1.0    ;   Sigw2=0.856   ;   Beta2=0.0828   ;   Gama2=0.44

!Part 2:
 Do I=1,NC
    WTN(1,I) = WTNP1(1,I)
    WTN(2,I) = WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
          
   !Part 4:  
	RKco=RKJ(NS)
          
   !Part 5:
    Call Kw_BC(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,NX,NY,DW,Mu,WB,WNP1,WTNP1,WTB)

   !Part 6:
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)

   !Part 7:
    Call KFi_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WTNP1,WTB,WB,DKX_C,DKY_C,DOmegX_C,DOmegY_C)
 
   !Part 8:
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,WNP1,WTNP1,WB,WTB,DKX_F,DKY_F,DOmegX_F,DOmegY_F)

   !part 9:
    Call KwSST_Func(Dim,NC,DW,WNP1,WTNP1,Mu,MR,DKX_C,DKY_C,DOmegX_C,DOmegY_C,Sigk1,Sigk2,&
                    Sigw1,Sigw2,Beta1,Beta2,Gama1,Gama2,Bstar,F11,F22,Sigk,Sigw,Beta,Gama)

   !Part 10:
	Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,WNP1,WTNP1,WB,WTB,Cont)
       
   !Part 11:
	Call KwSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DOmegX_F,DOmegY_F,MR,Sigk,Sigw,WNP1,WTNP1,Mu,Mut,Dift)
	               
   !Part 12:
	Call KwSST_Source(Dim,NC,MR,A,WNP1,WTNP1,Mut,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DOmegX_C,DOmegY_C,&
	                  Bstar,Sigw2,Beta,Gama,F11,St)

   !Part 13:
    Do J=1,NC

	   Co = RKco*DT(J)/A(J)

       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) - St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) - St(2,J) ) 

      !Part 14: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)

      !Part 15:
       Rho   = WNP1(1,J)
       K     = WTNP1(1,J)/Rho
       Omega = WTNP1(2,J)/Rho

       Vor= Dabs(DUY_C(J)-DVX_C(J))

      !Part 16:     
       Mut(J) = 0.31*Rho*K / max(0.31*Omega,Vor*F22(J)) /MR

    End Do !J
            
 END DO !NS
!*********************************************************************************************
 End
!###########################################################################################
