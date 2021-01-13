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
 Subroutine KFi_FaceGrad(Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB,DKX,DKY,DFiX,DFiY)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,Wnp1,WTNP1,Wb,WTB
 Intent(Out  )::DKX,DKY,DFiX,DFiY

 Integer::Dim,I,NP,NFW1,NFW2,NF,NF1,NF2,ME,NE,P1,P2
 Real(8)::K,Fi,DX,DY,AREA,KME,FiME,KNE,FiNE,DX1,DX2,DX3,DX4,DY1,DY2,DY3,DY4,DXC,DYC,K1,K2,Fi1,Fi2
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:Dim)::X,Y,XC,YC,DW,KP,FiP,DKX,DKY,DFiX,DFiY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NEC
!************************************************************************************************
!Part 1:
 DO I=1,NP
    KP(I)  = 0.0
    FiP(I) = 0.0
    NEC(I) = 0
 END DO

!Part 2:
 Do I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
	P1 = IDS(3,I)
    P2 = IDS(4,I)

   !Part 4:
    K  = 0.5 * ( WTNP1(1,ME)/WnP1(1,ME) + WTNP1(1,NE)/WnP1(1,NE) )
    Fi = 0.5 * ( WTNP1(2,ME)/WnP1(1,ME) + WTNP1(2,NE)/WnP1(1,NE) )
 
   !Part 5:
    KP(P1)  = KP(P1)  + K
    FiP(P1) = FiP(P1) + Fi
    
    KP(P2)  = KP(P2)  + K
    FiP(P2) = FiP(P2) + Fi

   !Part 6:
    NEC(P1) = NEC(P1) + 1
    NEC(P2) = NEC(P2) + 1

 End Do
 
!Part 7:
 Do I=NF2+1,NF
 
   !Part 8:
	P1 = IDS(3,I)
    P2 = IDS(4,I)
 
   !Part 9:
    K  = WTB(1,I)/WB(1,I)
    Fi = WTB(2,I)/WB(1,I)
 
   !Part 10:
    KP(P1)  = KP(P1)  + K 
    FiP(P1) = FiP(P1) + Fi 
    
    KP(P2)  = KP(P2)  + K
    FiP(P2) = FiP(P2) + Fi
      
    !Part 11:
    NEC(P1) = NEC(P1) + 1
    NEC(P2) = NEC(P2) + 1

 End Do

!Part 12:
 DO I=1,NP
    KP(I)  = KP(I) /NEC(I)
    FiP(I) = FiP(I)/NEC(I)
 END DO

!Part 13:
 DO I=NFW1+1,NFW2

    P1 = IDS(3,I)
    P2 = IDS(4,I)

    Fi = WTB(2,I)/WB(1,I)

    KP(P1)  = 0.0
    FiP(P1) = Fi
    
    KP(P2)  = 0.0
    FiP(P2) = Fi
 END DO
 
!Part 14:
 DO I=NF1+1,NF2
  
   !Part 15:        
    ME = IDS(1,I)        
	NE = IDS(2,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
  
   !Part 16:
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
  
   !Part 17:
    K1    = KP(P1)    
    K2    = KP(P2)  
    KME    = WTNP1(1,ME)/WnP1(1,ME)
    KNE    = WTNP1(1,NE)/WnP1(1,NE)

    Fi1  = FiP(P1)
    Fi2  = FiP(P2)
    FiME = WTNP1(2,ME)/WnP1(1,ME)
    FiNE = WTNP1(2,NE)/WnP1(1,NE)
    
    DKX(I) = 0.5*( (K1 + KNE)*DY1 + (KNE + K2)*DY2 + (KME + K2 )*DY3 + (K1 + KME)*DY4 )/AREA
    DKY(I) =-0.5*( (K1 + KNE)*DX1 + (KNE + K2)*DX2 + (KME + K2 )*DX3 + (K1 + KME)*DX4 )/AREA

    DFiX(I) = 0.5*( (Fi1 + FiNE)*DY1 + (FiNE + Fi2)*DY2 + (FiME + Fi2)*DY3 + (Fi1 + FiME)*DY4 )/AREA
    DFiY(I) =-0.5*( (Fi1 + FiNE)*DX1 + (FiNE + Fi2)*DX2 + (FiME + Fi2)*DX3 + (Fi1 + FiME)*DX4 )/AREA

 END DO
 
!Part 18:
 DO I=NF2+1,NF
  
   !Part 19:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
  
   !Part 20:
    DX1 = X (P2) - X (P1)
    DY1 = Y (P2) - Y (P1)
    DX2 = Xc(ME) - X (P1)
	DY2 = Yc(ME) - Y (P1)
    DX3 = Xc(ME) - X (P2)
    DY3 = Yc(ME) - Y (P2)
    DX4 = X (P1) - Xc(ME)
    DY4 = Y (P1) - Yc(ME)
    AREA = 0.5*Dabs( DX1*DY2 - DY1*DX2 )
  
   !Part 21:
    KME = WTNP1(1,ME)/WnP1(1,ME) 
    K1  = KP(P1)
    K2  = KP(P2)

    FiME = WTNP1(2,ME)/WnP1(1,ME)
    Fi1  = FiP(P1)
    Fi2  = FiP(P2)
   
    DKX(I)  = 0.5*( (KME  + K2 )*DY3  + (KME  + K1 )*DY4 + (K1  + K2 )*DY1 )/AREA 
    DKY(I)  =-0.5*( (KME  + K2 )*DX3  + (KME  + K1 )*DX4 + (K1  + K2 )*DX1 )/AREA
    DFiX(I) = 0.5*( (FiME + Fi2)*DY3  + (FiME + Fi1)*DY4 + (Fi1 + Fi2)*DY1 )/AREA  
    DFiY(I) =-0.5*( (FiME + Fi2)*DX3  + (FiME + Fi1)*DX4 + (Fi1 + Fi2)*DX1 )/AREA    
    
 END DO
!*********************************************************************************************
 End
!###########################################################################################

