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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
  Subroutine BC_Riemann_RoePreCond(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
 Implicit None                    
!*********************************************************************************************
 Intent (In   )::Dim,NFF1,NFF2,GM,U0,V0,P0,R0,C0,IDS,Wnp1,NX,NY,DA,P
 Intent (Out  )::Wb

 Integer::Dim,I,NFF1,NFF2,ME,P1,P2
 Real(8)::GM,GM1,C0,U0B,V0B,P0B,R0B,C0B,S0,DX,DY,DH,NXX,NYY,QN0,QT0,PE,RE,CE,&
          VE,QNE,QTE,SE,UE,QNN,C,QTT,SB,UB,VB,RB,PB,REB,U0,V0,P0,R0
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:5,1:Dim)::Wb
!*********************************************************************************************	
 GM1= GM-1.
 
!Part 1
 DO I=NFF1+1,NFF2

   !Part 2:
    ME  = IDS(1,I)
    P1  = IDS(3,I)
    P2  = IDS(4,I)
    NXX = NX(I)/DA(I)
    NYY = NY(I)/DA(I)

   !Part 3:
    U0B = U0
    V0B = V0
    P0B = P0
    R0B = R0
    C0B = SQRT(GM*P0B/R0B)

   !Part 4:
    RE = WNP1(1,ME)
    UE = WNP1(2,ME)/RE
    VE = WNP1(3,ME)/RE
    PE = P(ME)
    CE = SQRT(ABS(GM*PE/RE))
    
   !Part 5:
    QN0 = U0B*NXX+V0B*NYY
    QT0 =-U0B*NYY+V0B*NXX
    S0  = P0B/R0B**GM

   !Part 6:
    QNE = UE*NXX+VE*NYY
    QTE =-UE*NYY+VE*NXX
    SE  = PE/RE**GM
	      
    !Part 7:
     QNN = QN0
     C   = C0
     QTT = QT0
     SB  = S0
     PB=CE*CE*RE/GM

    !Part 8:
	 IF(QNN>0.0)then
      QNN = QNE
      C   = CE
      QTT = QTE
      SB  = SE
	  PB=C0*C0*R0/GM
	 Endif

   
   !Part 9:
    RB = (C*C/GM/SB)**(1.0/GM1)
    UB = QNN*NXX-QTT*NYY
    VB = QTT*NXX+QNN*NYY    
    REB= PB/GM1 + 0.5*RB*(UB*UB + VB*VB)

   !Part 10:
    Wb(1,I) = RB
    Wb(2,I) = RB*UB
    Wb(3,I) = RB*VB
    Wb(4,I) = REB
    Wb(5,I) = PB

 End do
!*********************************************************************************************
 End
!###########################################################################################