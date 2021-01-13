!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Turbulence Variables at Boundary Faces                    //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,H Kharinezhad Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F047F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Ke_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,WTB,MR,&
     DW,wnp1,mu,rokinf,roeinf)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,WB,WTNP1,IDS,MR,DW&
     ,WNP1,MU,rokinf,roeinf
 Intent(Out  )::WTB

 Real(8)::U,V,Q,MR,vn,rokinf,roeinf
 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME,P1,P2
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,DW,MU
 Real(8),Dimension(1:2,1:Dim)::WTB,WTNP1
 Real(8),Dimension(1:5,1:Dim)::WB
  Real(8),Dimension(1:4,1:Dim)::WNP1

!*******************************************************************************************

!Part 1:
 DO I=NFI1+1,NFI2
    
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    
 	WTB(1,I) = rokinf
    WTB(2,I) = roeinf

 END do

!Part 2:  
 DO I=NFO1+1,NFO2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 3:
 DO I=NFW1+1,NFW2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)  !!MR * WB(2,I)*(2.0*WTNP1(2,ME)*Mu(ME)/ (WNP1(1,ME)*DW(ME)*DW(ME)) ) !
    WTB(2,I) = WTNP1(2,ME)

 END do
  
!Part 4:
 DO I=NFS1+1,NFS2
    ME  = IDS(1,I)
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)
 END do
  
!Part 5:
 DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    
    Q  = U*NX(I)+V*NY(I)
        
    if (Q<=0.0) then
 	WTB(1,I) = rokinf
    WTB(2,I) = roeinf
    else
 	WTB(1,I) = WTNP1(1,ME)
    WTB(2,I) = WTNP1(2,ME)    
    end if

 END do
!*******************************************************************************************
 End
!###########################################################################################
