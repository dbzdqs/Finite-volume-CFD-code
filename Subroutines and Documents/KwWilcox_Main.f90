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
!// Developed by: H. Nazari, Mechanical Eng., Amirkabir University of Technology           //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KwWilcox_Main(Dim,NC,NP,DUY,INW,IDS,MR,NRKS,P,WTNP1,WNP1,WB,GM,NF,NF1,NF2,DA&
                     ,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,NX,NY,XC,YC,DW,A,Mu,DT,Mut)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NC,GM,MR,NRKS,P,WNP1,WB,DA,NFW1,NFW2,NFF1,NFF2&
                 ,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,X,Y,NX,NY,XC,YC,DW,A,Mu,DT,NP
 Intent (Out  )::WTNP1,Mut

Integer::Dim,I,J,NC,NF,NF1,NF2,NRKS,NP,NS,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,ME
Real(8)::RKco,CO,Omega,K,Rho,GM,ALFA,BETA,BETA_S,Sigk,Sigw,MR
Integer,Dimension(1:4,1:Dim)::IDS
Real(8),Dimension(1:2,1:Dim)::WTN,WTNP1,WTB,ConT,DifT,St
Real(8),Dimension(1:4,1:Dim)::WNP1
Real(8),Dimension(1:5,1:Dim)::WB
Real(8),Dimension(1:Dim)::Mut,P,DUX_C,DUY_C,DVX_C,DVY_C,DUY,DA
Real(8),Dimension(1:Dim)::X,Y,NX,NY,XC,YC,DW,A,Mu,DT,DKX,DKY,DWX,DWY
Integer,Dimension(1:Dim)::INW
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
	RKco=1.0/(NRKS-NS+1)
       
   !Part 5:
    Call Kw_BC(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,NX,NY,DW,Mu,WB,WNP1,WTNP1,WTB)

   !Part 6:  
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
            
   !Part 7:  
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX,DKY,DWX,DWY)

    Do I=NFS1+1,NFS2
       DKY(I)=0.0
	   DWY(I)=0.0
    End do

   !Part 8:
    Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,WNP1,WTNP1,WB,WTB,Cont)
  
   !Part 9:
    Call KFi_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX,DKY,DWX,DWY,Sigk,Sigw,MR,Mu,Mut,Dift)

   !Part 10:
    Call KwWilcox_Source(Dim,NC,MR,ALFA,BETA,BETA_S,Mut,A,WTNP1,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C,St)
      
   !Part 11:
    Do J=1,NC

	   Co = RKco*DT(J)/A(J)
        
       WTNP1(1,J) = WTN(1,J)  - Co*( ConT(1,J) + DifT(1,J) - St(1,J) ) 
       WTNP1(2,J) = WTN(2,J)  - Co*( ConT(2,J) + DifT(2,J) - St(2,J) )
	    
      !Part 12: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)

      !Part 13:
       Rho       = WNP1(1,J)
       K         = WTNP1(1,J)/Rho
       Omega     = WTNP1(2,J)/Rho
         
      !Part 15:  
       Mut(J) =(1.0/MR)*Rho*K / (Omega+1.e-20)

   End Do !J

 End Do !NS
  
!*********************************************************************************************
 End
!###########################################################################################