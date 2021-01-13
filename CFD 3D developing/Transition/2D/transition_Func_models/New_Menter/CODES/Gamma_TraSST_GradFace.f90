!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// calculation of Face Gradients for transition model                                   //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!**********************************************************************************************
 Subroutine Gamma_TraSST_GradFace(Dim,NC,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,DGX,DGY,DW,MR,&
                                  beta1,Wnp1,Wntp1,Wb,Wbt,Mu,Mut)
 Implicit None
!**********************************************************************************************
 Integer::Dim,I,NP,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE,P1,P2
 Real(8)::K,Omeg,G,DX,DY,AREA,&
          KME,OmegME,KNE,OmegNE,GME,GNE,&
		  DX1,DX2,DX3,DX4,DY1,DY2,DY3,DY4,DXC,DYC,K1,K2,Omeg1,Omeg2,G1,G2,beta1,MR
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,DW,KP,OmegP,GP,Mu,Mut,DKX,DKY,DOmegX,DOmegY,DGX,DGY
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Wbt
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NEC
!***************************************************************************************************
 !Part 1:
 DO I=1,NF
    
    DGX(I) =  0.0
    DGY(I) =  0.0
    
 END DO
 
 !Part 2:
 DO I=1,NP
    
    GP(I) = 0.0
    
    NEC(I)= 0
 END DO

 
!Part 3:
 Do I=NF1+1,NF2
   !Part 4:
    ME = IDS(1,I)
    NE = IDS(2,I)
	P1 = IDS(3,I)
    P2 = IDS(4,I)

   !Part 5:
    
    G    = 0.5 * ( WntP1(3,ME)/WnP1(1,ME) + WntP1(3,NE)/WnP1(1,NE) )
 
   !Part 6:
    GP(P1) = GP(P1) + G
    
    GP(P2) = GP(P2) + G
   
    !Part 7:
    NEC(P1) = NEC(P1) + 1
    NEC(P2) = NEC(P2) + 1

 End Do
 
 
 !Part 8:
 Do I=NFW2+1,NF
 
   !Part 9:
	P1 = IDS(3,I)
    P2 = IDS(4,I)
 
   !Part 10:
    G    = Wbt(3,I)/Wb(1,I)
    
 
   !Part 11:
    GP(P1) = GP(P1) + G

    GP(P2) = GP(P2) + G
 
    !Part 12:
    NEC(P1) = NEC(P1) + 1
    NEC(P2) = NEC(P2) + 1

 End Do

!Part 13:
 DO I=1,NP
    KP(I) = KP(I)/NEC(I)
    OmegP(I) = OmegP(I)/NEC(I)
    
    GP(I) = GP(I)/NEC(I)
    
 END DO

  !Part 14:
 DO I=NFW1+1,NFW2
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
    
    GP(P1) = WntP1(3,ME)/WnP1(1,ME) 
    
    GP(P2) = WntP1(3,ME)/WnP1(1,ME) 
    
 END DO 
 
!Part 15:
 DO I=NF1+1,NF2
  
   !Part 16:        
    ME = IDS(1,I)        
	NE = IDS(2,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
  
   !Part 17:
    DX1 = XC(NE) - X (P1)
    DX2 = X (P2) - XC(NE)
    DX3 = XC(ME) - X (P2)
    DX4 = X (P1) - XC(ME)
    DY1 = YC(NE) - Y (P1)
    DY2 = Y (P2) - YC(NE)
    DY3 = YC(ME) - Y (P2)
    DY4 = Y (P1) - YC(ME)

    DX   = X(P2)  - X(P1)
    DY   = Y(P2)  - Y(P1)
    DXC  = XC(ME) - XC(NE)
    DYC  = YC(ME) - YC(NE)
    AREA = ABS(0.5*(DX*DYC - DY*DXC))
  
  
   !Part 18:
    G1=GP(P1)
    
    G2=GP(P2)
    
    GME = WntP1(3,ME)/WnP1(1,ME)
  
    GNE = WntP1(3,NE)/WnP1(1,NE)
        
    DGX(I) = 0.5*( (G1 + GNE)*DY1 + (GNE + G2)*DY2 + (GME + G2 )*DY3 + (G1  + GME)*DY4 )/AREA
    DGY(I) =-0.5*( (G1 + GNE)*DX1 + (GNE + G2)*DX2 + (GME + G2 )*DX3 + (G1  + GME)*DX4 )/AREA
    

 END DO

!Part 19:
 DO I=NFW1+1,NFW2
  
   !Part 20:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
  
   !Part 21:
    DX = X (P2) - X (P1)
    DY = Y (P2) - Y (P1) 
    DX1 = X (P2)-X(P1)
	DY1 = Y (P2)-Y(P1)
    DX2 = Xc(ME)-X(P1)
	DY2 = Yc(ME)-Y(P1)
    DX3 = XC(ME) - X (P2)
    DY3 = YC(ME) - Y (P2)
    DX4 = X (P1) - XC(ME)
    DY4 = Y (P1) - YC(ME)
    AREA = 0.5*Dabs( DX1*DY2 - DY1*DX2 )
  
   !Part 22:  
    
    GME = WntP1(3,ME)/WnP1(1,ME)
    G1=GP(P1)
    G2=GP(P2)
        
    DGX(I) = 0.5*( (GME + G2)*DY3  + (GME + G1)*DY4 + (G1 + G2)*DY1 )/AREA  
    DGY(I) = -0.5*( (GME + G2)*DX3  + (GME + G1)*DX4 + (G1 + G2)*DX1 )/AREA    
    
              
 END DO
 
    !Part 23:
  DO I=NFW2+1,NF
  
   !Part 24:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
  
   !Part 25:
    DX = X (P2) - X (P1)
    DY = Y (P2) - Y (P1)
    DX1 = X (P2)-X(P1)
    DY1 = Y (P2)-Y(P1)
    DX2 = Xc(ME)-X(P1)
	DY2 = Yc(ME)-Y(P1)
    DX3 = XC(ME) - X (P2)
    DY3 = YC(ME) - Y (P2)
    DX4 = X (P1) - XC(ME)
    DY4 = Y (P1) - YC(ME)
    AREA = 0.5*Dabs( DX1*DY2 - DY1*DX2 )
  

   !Part 26:
    
    GME = WntP1(3,ME)/WnP1(1,ME) 
    G1 = GP(P1)
    G2 = GP(P2)
       
    DGX(I) = 0.5*( (GME + G2)*DY3  + (GME + G1)*DY4 + (G1 + G2)*DY1 )/AREA  
    DGY(I) =-0.5*( (GME + G2)*DX3  + (GME + G1)*DX4 + (G1 + G2)*DX1 )/AREA     
           
  END DO
 

!**********************************************************************************************
 End
!##############################################################################################

