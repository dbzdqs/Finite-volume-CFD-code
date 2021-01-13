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
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KeYang_Main(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                         NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,&
                NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,I,II,J,NC,NS,NP,NRKS,ME,NF,NFW1,NFW2,NF1,NF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::RKco,K,Epsilon,Rho,MR,Co,tauwall,ustar,yplus,Yn,Sigk,Sige,Ce1,Ce2,Cmu
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,WTN,Cont,Dift,St
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,XC,YC,DW,A,Mu,Mut,DT,DUY
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DEPSX_C,DEPSY_C,DD2UY
 Real(8),Dimension(1:Dim)::DKX_F,DKY_F,DEPSX_F,DEPSY_F
 Real(8),Dimension(1:Dim)::Lk,Le,f1,f2,fmu
!*********************************************************************************************
!Part 1:
 sigk=1.0d0
 sige=1.0/1.3
 Ce1=1.35d0
 Ce2=1.8d0
 Cmu=0.09d0 

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
    Call KeYang_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,Wnp1,WTNP1,IDS,MR,Mu,WTB)
         
   !Part 6:  
    Call Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
   
    Call KFi_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WTNP1,WTB,WB,DKX_C,DKY_C,DEPSX_C,DEPSY_C)

    Call Grad2AtCell(Dim,NC,NF,NFW1,NF1,NF2,IDS,A,NX,NY,DUY_C,DD2UY)
    
   !Part 7:  
    Call KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX_F,DKY_F,DEPSX_F,DEPSY_F)
    Do I=NFS1+1,NFS2
       DKY_F(I)=0.0
	   DEPSY_F(I)=0.0
    End do
      
   !Part 8:
    Call KeYang_Funcs(Dim,NC,IDS,X,Y,Xc,Yc,NX,NY,A,DW,INW,MR,Wnp1,WTNP1,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DEPSX_C,DEPSY_C,&
            DD2UY,Mu,Mut,f1,f2,fmu,Lk,Le)
         
   !Part 9:
	Call KFi_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,WNP1,WTNP1,WB,WTB,Cont)

   !Part 10:
	Call KFi_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,MR,Mu,Mut,Dift)
         
   !Part 11:
    Call KeYang_Source(Dim,IDS,NC,A,MR,Wnp1,WTNP1,Mu,Mut,DUX_C,DUY_C,DVX_C,DVY_C,DKX_C,DKY_C,DEPSX_C,DEPSY_C,Ce1,Ce2,f1,f2,Lk,Le,St)
         
   !Part 12:
    Do J=1,NC

	   Co = RKco*DT(J)/A(J)
         
       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) + St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) + St(2,J) ) 

      !Part 13: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)
	
      !Part 14:
       Rho     = Wnp1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho
		       
      !Part 15:  
       Mut(J) = (1.0/MR)*Cmu*Fmu(J)*Rho*k*k / Epsilon 
              
      End Do !J
      
 End Do !NS
!*********************************************************************************************
 End
!###########################################################################################
