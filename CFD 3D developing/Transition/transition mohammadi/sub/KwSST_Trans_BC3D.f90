!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Turbulence and transition Variables at Boundary Faces     //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar,M.A.Zoljanahi Iran, Tehran, OpenFlows@chmail.ir              //!
!// Doc ID: MC2F057F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_BC3D(Dim,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,&
                             MR,NX,NY,NZ,DW,Mu,WB,WNP1,WTNP1,Kinf,Oinf,Ginf,WTB)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NX,NY,NZ,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,IDS,MR,DW,&
                 Mu,WB,WNP1,WTNP1,Kinf,Oinf,Ginf
 Intent (Out  )::WTB

 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME
 Real(8)::MR,U,V,W,Q,Kinf,Oinf,Ginf
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Mu,DW,NX,NY,NZ
 Real(8),Dimension(1:3,1:Dim)::WTB,WTNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:Dim)::WNP1
!*******************************************************************************************

!Part 1:
 DO I=NFI1+1,NFI2
 	WTB(1,I) = Kinf
    WTB(2,I) = Oinf 
    WTB(3,I) = Ginf
 END do
  
!Part 2:  
 DO I=NFO1+1,NFO2
    ME       = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
    WTB(3,I) = WTNP1(3,ME)
 END do
  
!Part 3:
 DO I=NFW1+1,NFW2
    ME       = IDS(1,I)
 	WTB(1,I) = 0.0
    WTB(2,I) = MR * WB(1,I)*(60.0*Mu(ME)/ (WNP1(1,ME)*0.075*DW(ME)*DW(ME)) ) 
    WTB(3,I) = WTNP1(3,ME)
 END do
 
!Part 4:
 DO I=NFS1+1,NFS2
    ME       = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
    WTB(3,I) = WTNP1(3,ME)
 END do
  
!Part 5:
 DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)
    
    Q  = U*NX(I)+V*NY(I)+W*NZ(I)
        
    IF( Q<=0. )Then  
 	 WTB(1,I) = Kinf
     WTB(2,I) = Oinf
     WTB(3,I) = Ginf
    Else
 	 WTB(1,I) = WTNP1(1,ME)
     WTB(2,I) = WTNP1(2,ME)
     WTB(3,I) = WTNP1(3,ME)     
    End IF
 END Do
!*******************************************************************************************
 End
!###########################################################################################

