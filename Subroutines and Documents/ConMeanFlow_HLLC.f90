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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: S. Kavoosi, Mechanical Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_HLLC(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L,E_L
 Real(8)::R_R,U_R,V_R,H_R,P_R,E_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 
 
 Real(8):: R_BAR,A_BAR,P_BAR,P_PVRS,P_STAR,Q_L,Q_R,S_L,S_R,S_STAR,UL_COEF,UR_COEF
 Real(8):: U1_SL,U2_SL,U3_SL,U4_SL,U1_SR,U2_SR,U3_SR,U4_SR
 Real(8):: F1_SL,F2_SL,F3_SL,F4_SL,F1_SR,F2_SR,F3_SR,F4_SR
 Real(8):: F1_HLLC,F2_HLLC,F3_HLLC,F4_HLLC
 
 Real(8),DIMENSION(1:4)::DEL_W
 Real(8):: G1_HLLC,G2_HLLC,G3_HLLC,G4_HLLC
 Real(8)::G1_L,G2_L,G3_L,G4_L
 Real(8)::G1_R,G2_R,G3_R,G4_R
 
 Real(8)::P_SL,P_SR,ALP_L,BET_L,OMG_L,ALP_R,BET_R,OMG_R,Z_TR,P_S
 integer::ISTR
 
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I) + V*NY(I)
    Pm = WB(5,I)

   !Part 6:
    F1 = Q * Wb(1,I) 
    F2 = Q * Wb(2,I) + Pm*NX(I)
    F3 = Q * Wb(3,I) + Pm*NY(I)
    F4 = Q *(Wb(4,I)+Pm)

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    L = IDS(1,I)
    R = IDS(2,I)

   !Part 10:
    DAA = DA(I)
	NXX = NX(I)/DAA
    NYY = NY(I)/DAA

   !Part 11:
    R_L = WNP1(1,L)
    U_L = WNP1(2,L)/WNP1(1,L)
    V_L = WNP1(3,L)/WNP1(1,L)
    P_L = P(L)
    E_L = WNP1(4,L)/WNP1(1,L) 
    
    R_R = WNP1(1,R)
    U_R = WNP1(2,R)/WNP1(1,R)
    V_R = WNP1(3,R)/WNP1(1,R)
    P_R = P(R)
    E_R = WNP1(4,R)/WNP1(1,R) 
    
   !Part 12: 
    H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
    H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
    
    A_L  = DSQRT(GM*P(L)/WNP1(1,L))
	A_R  = DSQRT(GM*P(R)/WNP1(1,R))
    
   !Part 13: 
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
   !Part 14:
    ISTR = 5
   
   !Part 15: 
    IF (ISTR==1) THEN
        
        S_L = UC_L-A_L
        S_R = UC_R+A_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
   !Part 16:     
    ELSE IF (ISTR==2) THEN
        
        S_L = DMIN1(UC_L-A_L,UC_R-A_R)
        S_R = DMAX1(UC_L+A_L,UC_R+A_R)
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
   
   !Part 17:     
    ELSE IF (ISTR==3) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = UC_ROE-A_ROE
        S_R = UC_ROE+A_ROE
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
   !Part 18:     
    ELSE IF (ISTR==4) THEN
        
        R_RT = DSQRT(R_R/R_L)
        R_ROE = DSQRT(R_L*R_R)
        U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
        V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
        H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
        A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
        
        UC_ROE = U_ROE*NXX+V_ROE*NYY
        
        S_L = DMIN1(UC_L-A_L,UC_ROE-A_ROE)
        S_R = DMAX1(UC_R+A_R,UC_ROE+A_ROE)
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
   !Part 19:     
    ELSE IF (ISTR==5) THEN
        
        R_BAR = 0.5*(R_L+R_R)                                                                        
        A_BAR = 0.5*(A_L+A_R)
        P_BAR = 0.5*(P_L+P_R)
        
        P_PVRS = P_BAR-0.5*(UC_R-UC_L)*R_BAR*A_BAR
        P_STAR = DMAX1(0.D0,P_PVRS)
        
        IF (P_STAR<=P_L) THEN
            
            Q_L = 1.0
            
        ELSE
            
            Q_L = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_L)-1.0))
            
        END IF
        
        IF (P_STAR<=P_R) THEN
            
            Q_R = 1.0
            
        ELSE
            
            Q_R = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_R)-1.0))
            
        END IF
        
        S_L = UC_L-A_L*Q_L
        S_R = UC_R+A_R*Q_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
    !Part 20:     
    ELSE IF (ISTR==6) THEN
        
        Z_TR = (GM-1.0)/(2.0*GM)
        
        P_STAR = ((A_L+A_R-0.5*(GM-1.0)*(U_R-U_L))/((A_L/(P_L**Z_TR))+(A_R/(P_R**Z_TR))))**(1.0/Z_TR)
       
        IF (P_STAR<=P_L) THEN
            
            Q_L = 1.0
            
        ELSE
            
            Q_L = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_L)-1.0))
            
        END IF
        
        IF (P_STAR<=P_R) THEN
            
            Q_R = 1.0
            
        ELSE
            
            Q_R = DSQRT(1.0+((GM+1.0)/(2*GM))*((P_STAR/P_R)-1.0))
            
        END IF
        
        S_L = UC_L-A_L*Q_L
        S_R = UC_R+A_R*Q_R
        
        S_STAR = (P_R-P_L+R_L*UC_L*(S_L-UC_L)-R_R*UC_R*(S_R-UC_R))/(R_L*(S_L-UC_L)-R_R*(S_R-UC_R))
        
    END IF
     
    !Part 21:     
     P_SL = P_L+R_L*(UC_L-S_L)*(UC_L-S_STAR)
     P_SR = P_R+R_R*(UC_R-S_R)*(UC_R-S_STAR)
     
     P_S = 0.5*(P_SL+P_SR)
     
    !Part 22:
     ALP_L = (S_L-UC_L)/(S_L-S_STAR)
     BET_L = ALP_L*(S_STAR-UC_L)
    
     ALP_R = (S_R-UC_R)/(S_R-S_STAR)
     BET_R = ALP_R*(S_STAR-UC_R)
     
    !Part 23: 
     U1_SL = ALP_L*(WNP1(1,L)) + (0.0)
     U2_SL = ALP_L*(WNP1(2,L)) + (R_L*BET_L*NXX)
     U3_SL = ALP_L*(WNP1(3,L)) + (R_L*BET_L*NYY)
     U4_SL = ALP_L*(WNP1(4,L)) + ((P_S*S_STAR-P_L*UC_L)/(S_L-S_STAR))
     
     U1_SR = ALP_R*(WNP1(1,R)) + (0.0)
     U2_SR = ALP_R*(WNP1(2,R)) + (R_R*BET_R*NXX)
     U3_SR = ALP_R*(WNP1(3,R)) + (R_R*BET_R*NYY)
     U4_SR = ALP_R*(WNP1(4,R)) + ((P_S*S_STAR-P_R*UC_R)/(S_R-S_STAR))
     
    !Part 24: 
     F1_L = (UC_L*R_L               )
     F2_L = (UC_L*WNP1(2,L)+P(L)*NXX)
     F3_L = (UC_L*WNP1(3,L)+P(L)*NYY)
     F4_L = (UC_L*(WNP1(4,L)+P(L))  )
    
     F1_R = (UC_R*R_R               )
     F2_R = (UC_R*WNP1(2,R)+P(R)*NXX)
     F3_R = (UC_R*WNP1(3,R)+P(R)*NYY)
     F4_R = (UC_R*(WNP1(4,R)+P(R))  )
     
    !Part 25: 
     IF (S_L>=0.0) THEN
         
         F1 = F1_L
         F2 = F2_L
         F3 = F3_L
         F4 = F4_L
         
    !Part 26:     
     ELSE IF (S_L<=0.0 .AND. S_STAR>=0.0) THEN
         
         F1 = F1_L+S_L*(U1_SL-WNP1(1,L))
         F2 = F2_L+S_L*(U2_SL-WNP1(2,L))
         F3 = F3_L+S_L*(U3_SL-WNP1(3,L))
         F4 = F4_L+S_L*(U4_SL-WNP1(4,L))
         
    !Part 27:     
     ELSE IF (S_STAR<=0.0 .AND. S_R>=0.0) THEN
         
         F1 = F1_R+S_R*(U1_SR-WNP1(1,R))
         F2 = F2_R+S_R*(U2_SR-WNP1(2,R))
         F3 = F3_R+S_R*(U3_SR-WNP1(3,R))
         F4 = F4_R+S_R*(U4_SR-WNP1(4,R))
         
    !Part 28:     
     ELSE IF (S_R<=0) THEN
         
         F1 = F1_R
         F2 = F2_R
         F3 = F3_R
         F4 = F4_R
         
     END IF
    
   !Part 29: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 30:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End 
!###########################################################################################
