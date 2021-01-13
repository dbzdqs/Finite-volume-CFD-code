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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BC_RiemannForShockTube(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,WNP1,P,WB,XC)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NFF1,NFF2,GM,IDS,Wnp1,NX,NY,DA,P,XC
 Intent (Out  )::Wb
 

 Integer::Dim,I,NFF1,NFF2,ME,P1,P2
 Real(8)::GM,GM1,U,V,CC,MLoc,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,STH,CTH,QN0,QT0,RI0,PE,RE,CE,&
          VE,QNE,QTE,RIE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P,XC
 Real(8),Dimension(1:5,1:Dim)::Wb
!*********************************************************************************************	
  
 GM1= GM-1.

 DO I=NFF1+1,NFF2

    ME  = IDS(1,I)
    P1  = IDS(3,I)
    P2  = IDS(4,I)

    STH = NX(I)/DA(I)
    CTH = NY(I)/DA(I)

    U = WNP1(2,ME)/WNP1(1,ME)
    V = WNP1(3,ME)/WNP1(1,ME)

    CC = GM*P(ME)/WNP1(1,ME)

    MLoc = SQRT((U*U+V*V)/CC)
    
    
    IF(XC(me)<0.3)Then                       !!!OOOOOOOOOOOKKKKKKKKKKKKKKKK!!!
     R0 = 1.0 
     U0 = 0.75
     V0 = 0.0
     P0 = 1.0
     C0 = GM*P0/R0
    Else
     R0 = 0.125 
     U0 = 0.0
     V0 = 0.0
     P0 = 0.1
     C0 = GM*P0/R0
    Endif
    

    U0B = U0
    V0B = V0
    P0B = P0
    R0B = R0
    C0B = SQRT(GM*P0B/R0B)

    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
    CE = SQRT(ABS(GM*PE/RE))

    QN0 = U0B*STH+V0B*CTH
    QT0 =-U0B*CTH+V0B*STH
    RI0 = QN0 - 2.*C0B/GM1
    S0  = P0B/R0B**GM

    QNE = UE*STH+VE*CTH
    QTE =-UE*CTH+VE*STH
    RIE = QNE+2.*CE/GM1
    SE  = PE/RE**GM

    IF(MLoc<=1.)then

     QNN = (RIE+RI0)/2.0
     C   = 0.25*GM1*(RIE-RI0)
     QTT = QT0
     SB  = S0
     IF(QNN>0.0)then
	  QTT = QTE
      SB  = SE
	 Endif

    Elseif(MLoc>1.)then

     QNN = QN0
     C   = C0
     QTT = QT0
     SB  = S0
     
	 IF(QNN>0.0)then
      QNN = QNE
      C   = CE
      QTT = QTE
      SB  = SE
	 Endif

	Endif

    RB = (C*C/GM/SB)**(1.0/GM1)
    UB = QNN*STH-QTT*CTH
    VB = QTT*STH+QNN*CTH
    PB = C*C*RB/GM

    REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)

    Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

 End do

!*********************************************************************************************
 End
!###########################################################################################
