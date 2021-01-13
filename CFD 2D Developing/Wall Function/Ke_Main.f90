!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate Turbulence Viscosity and Fluctuation Term of Turbulent Flow//!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,H Kharinezhad Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC5F050F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Ke_Mainw(Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
                         NFF1,NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,WTNP1,Mut,IWF,p ,U0,V0,Mut0,RKJ)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,&
                NFF2,NP,IDS,INW,XC,YC,DW,DT,A,MR,NRKS,Wb,WNP1,Mu,DUY,p,U0,V0,Mut0
 Intent(Out  )::WTNP1,Mut
 Intent(inOut  ):: IWF
 

 Integer::Dim,I,II,J,NC,NS,NP,NRKS,ME,NF,NFW1,NFW2,NF1,NF2,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2
 Real(8)::RKco,K,Epsilon,Rho,MR,Co,tauwall,ustar,yplus,Yn,fmu,Sigk,Sige,Ce1,Ce2,Cmu, rm,rokinf,roeinf,U0,V0,Tu,Mut0
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,WTN,Cont,Dift,St
 Real(8),Dimension(1:3,1:Dim)::IWF
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,XC,YC,DW,A,Mu,Mut,DT,DUY,DUX_C,DUY_C,DVX_C,DVY_C,&
                           DKX_F,DKY_F,DEPSX_F,DEPSY_F,p
 Real(8)::c1,c2,kapa,kplus,eplus,y11(1:100),miu,y10(1:100),y12(1:100),y13(1:100)
 
 Real(8),Dimension(1:5)::RKJ
!*******************************************************************************************
Tu=.01
!Part 1:
 Sigk=1.0d0
 Sige=1.3
 Ce1=1.35d0
 Ce2=1.8d0
 Cmu=0.09d0
 rokinf=0.000025*(U0*U0+V0*V0) 
 roeinf=cmu*rokinf*rokinf/(Mut0*Tu) /MR
 
 
!Part 2:
 Do I=1,NC
    WTN(1,I) =WTNP1(1,I)
    WTN(2,I) =WTNP1(2,I)
 End Do

!Part 3:
 Do NS=1,NRKS
      
   !Part 4:    
	RKco=RKJ(NS)
       
   !Part 5:
    Call Ke_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB,MR,DW,wnp1,mu,rokinf,roeinf)

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
    Call Ke_Source(Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DUX_C,DUY_C,DVX_C,DVY_C,St)
    
        
    !Part 11:
    Call wallfunc_swfke(Dim,NC,NX,NY,NFW1,NFW2,IDS,DW,A,INW,MR,WNP1,WTNP1,Mu,Mut,iwf ,x,y,St)
    
       !Part 15:      
    Do J=NFW1+1,NFW2
        
       ME      = IDS(1,J)
       Rho     = WNP1(1,ME)
       K       = WTNP1(1,ME)/Rho
       Epsilon = WTNP1(2,ME)/Rho
  
       Mut(ME) = (1.0/MR)*Cmu*Rho*k*k / Epsilon 
     
    End Do !J
    
   !Part 12:
    Do J=1,NC

	   Co = RKco*DT(J)/A(J)
         
       WTNP1(1,J) = WTN(1,J) - Co*( Cont(1,J) + Dift(1,J) + St(1,J) ) 
       WTNP1(2,J) = WTN(2,J) - Co*( Cont(2,J) + Dift(2,J) + St(2,J) ) 

      !Part 13: 
       if(WTNP1(1,J)<0.0 ) WTNP1(1,J)=WTN(1,J)
       if(WTNP1(2,J)<0.0 ) WTNP1(2,J)=WTN(2,J)
	
      !Part 14:
       Rho     = WNP1(1,J)
       K       = WTNP1(1,J)/Rho
       Epsilon = WTNP1(2,J)/Rho

      !Part 15:  
       Mut(J) = (1.0/MR)*Cmu*Rho*k*k / Epsilon 
              
    End Do !J

    
End Do !NS
    
 
!*******************************************************************************************
 End
!###########################################################################################
