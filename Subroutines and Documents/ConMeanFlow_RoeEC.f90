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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_RoeEC(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
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
 
 Real(8)::R_L,U_L,V_L,H_L,P_L
 Real(8)::R_R,U_R,V_R,H_R,P_R

 Real(8)::UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 
 Real(8)::Z1_BAR,Z2_BAR,Z3_BAR,Z4_BAR,Z1_LN,Z4_LN
 Real(8)::R_HAT,U_HAT,V_HAT,P1_HAT,P2_HAT,A_HAT,H_HAT,UC_HAT
 Real(8)::F1_CNV,F2_CNV,F3_CNV,F4_CNV
 Real(8),DIMENSION(1:4,1:4)::R_REC,RT_REC,JCB_REC
 Real(8),DIMENSION(1:4)::EV_REC,S_REC,DV_REC,DSP_REC
 INTEGER::IREC,JREC,KREC
 
 Real(8)::LOG_AVE
 Real(8)::C1,C2
 
 INTEGER::IES
 Real(8)::DM_MAX,BET_EC2,ALF_MIN,ALF_MAX,ALF_EC2,ALF_EC1
 
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
    
    R_R = WNP1(1,R)
    U_R = WNP1(2,R)/WNP1(1,R)
    V_R = WNP1(3,R)/WNP1(1,R)
    P_R = P(R)
    
   !Part 12: 
    H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
    H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
    
    A_L  = DSQRT(GM*P(L)/WNP1(1,L))
	A_R  = DSQRT(GM*P(R)/WNP1(1,R))
    
    
   !Part 13:
    Z1_BAR = 0.5*(DSQRT(R_L/P_L)+DSQRT(R_R/P_R))
    Z2_BAR = 0.5*(DSQRT(R_L/P_L)*U_L+DSQRT(R_R/P_R)*U_R)
    Z3_BAR = 0.5*(DSQRT(R_L/P_L)*V_L+DSQRT(R_R/P_R)*V_R)              
    Z4_BAR = 0.5*(DSQRT(R_L*P_L)+DSQRT(R_R*P_R))                  
   
   !Part 14: 
    C1 = DSQRT(R_L/P_L)
    C2 = DSQRT(R_R/P_R)
    Z1_LN = LOG_AVE(C1,C2)
    
    C1 = DSQRT(R_L*P_L)
    C2 = DSQRT(R_R*P_R)
    Z4_LN = LOG_AVE(C1,C2)
    
   !Part 15: 
    R_HAT = Z1_BAR*Z4_LN
    U_HAT = Z2_BAR/Z1_BAR
    V_HAT = Z3_BAR/Z1_BAR
    P1_HAT = Z4_BAR/Z1_BAR
    P2_HAT = ((GM+1)*(Z4_LN))/(2.0*GM*Z1_LN)+((GM-1)*(Z4_BAR))/(2.0*GM*Z1_BAR)
    
    A_HAT = DSQRT((GM*P2_HAT)/R_HAT)
    H_HAT = (A_HAT**2)/(GM-1.0)+0.5*(U_HAT**2+V_HAT**2)
    
    UC_HAT = U_HAT*NXX+V_HAT*NYY

   !Part 16:  
    F1_CNV = UC_HAT*R_HAT
    F2_CNV = UC_HAT*R_HAT*U_HAT+P1_HAT*NXX
    F3_CNV = UC_HAT*R_HAT*V_HAT+P1_HAT*NYY
    F4_CNV = UC_HAT*R_HAT*H_HAT
    
   !Part 17:  
    R_REC(1,1) = 1.0
    R_REC(1,2) = 1.0
    R_REC(1,3) = 0.0
    R_REC(1,4) = 1.0
    
    R_REC(2,1) = U_HAT-A_HAT*NXX
    R_REC(2,2) = U_HAT
    R_REC(2,3) = NYY
    R_REC(2,4) = U_HAT+A_HAT*NXX
    
    R_REC(3,1) = V_HAT-A_HAT*NYY
    R_REC(3,2) = V_HAT
    R_REC(3,3) =-NXX
    R_REC(3,4) = V_HAT+A_HAT*NYY
    
    R_REC(4,1) = H_HAT-A_HAT*UC_HAT
    R_REC(4,2) = 0.5*(U_HAT**2+V_HAT**2)
    R_REC(4,3) = U_HAT*NYY-V_HAT*NXX
    R_REC(4,4) = H_HAT+A_HAT*UC_HAT
    
   !Part 18:  
    IES=2
    IF (IES==0) THEN
        
        !Part 19: 
        !!ENTROPY STABLE SCHEME, ROE-ES
        EV_REC(1) = DABS(UC_HAT-A_HAT)
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = DABS(UC_HAT+A_HAT)
        
    ELSE IF (IES==1) THEN
        
        !Part 20: 
        !!ENTROPY CONSISTENT#1 SCHEME, ROE-EC1
        UC_L = U_L*NXX+V_L*NYY
        UC_R = U_R*NXX+V_R*NYY
        
        !Part 21:
        ALF_EC1 = 1.0/6.0       
        
        EV_REC(1) = DABS(UC_HAT-A_HAT)+(ALF_EC1)*DABS( (UC_R-A_R)-(UC_L-A_L) )
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = DABS(UC_HAT+A_HAT)+(ALF_EC1)*DABS( (UC_R+A_R)-(UC_L+A_L) )
        
    ELSE IF (IES==2) THEN
        
        !Part 22: 
        !!ENTROPY CONSISTENT#2 SCHEME, ROE-EC2
        ALF_MIN = 1.0/6.0
        ALF_MAX = 2.0
        DM_MAX = 0.5
        BET_EC2 = 1.0/6.0
        
        !Part 23:
        UC_L = U_L*NXX+V_L*NYY
        UC_R = U_R*NXX+V_R*NYY
        
        !Part 24:
        M_L = DSQRT( (U_L**2+V_L**2)/(GM*P_L/R_L) )
        M_R = DSQRT( (U_R**2+V_R**2)/(GM*P_R/R_R) ) 
        
        !Part 25: 
        ALF_EC2 = (ALF_MAX-ALF_MIN)*( DMAX1(0.0,DSIGN(DM_MAX,(M_R-M_L))) ) + ALF_MIN               
        
        EV_REC(1) = (1.0+BET_EC2)*DABS(UC_HAT-A_HAT)+(ALF_EC2)*DABS( (UC_R-A_R)-(UC_L-A_L) )
        EV_REC(2) = DABS(UC_HAT)
        EV_REC(3) = DABS(UC_HAT)  
        EV_REC(4) = (1.0+BET_EC2)*DABS(UC_HAT+A_HAT)+(ALF_EC2)*DABS( (UC_R+A_R)-(UC_L+A_L) )
        
        
    END IF
    
    !Part 26:
    S_REC(1) = (R_HAT)/(2.0*GM)
    S_REC(2) = ((GM-1)*R_HAT)/GM
    S_REC(3) = P1_HAT                              
    S_REC(4) = (R_HAT)/(2.0*GM)
    
    !Part 27:
    DV_REC(1) = ( (GM-LOG(P_R/(R_R**GM)))/(GM-1.0)-0.5*(R_R/P_R)*(U_R**2+V_R**2) ) - ( (GM-LOG(P_L/(R_L**GM)))/(GM-1.0)-0.5*(R_L/P_L)*(U_L**2+V_L**2) )
    DV_REC(2) = ( (R_R/P_R)*U_R ) - ( (R_L/P_L)*U_L )
    DV_REC(3) = ( (R_R/P_R)*V_R ) - ( (R_L/P_L)*V_L )
    DV_REC(4) = -(R_R/P_R) + (R_L/P_L)
    
    !Part 28:
    RT_REC(1,1) = 1.0                                                 
    RT_REC(1,2) = U_HAT-A_HAT*NXX
    RT_REC(1,3) = V_HAT-A_HAT*NYY
    RT_REC(1,4) = H_HAT-A_HAT*UC_HAT
    
    RT_REC(2,1) = 1.0
    RT_REC(2,2) = U_HAT
    RT_REC(2,3) = V_HAT
    RT_REC(2,4) = 0.5*(U_HAT**2+V_HAT**2)
    
    RT_REC(3,1) = 0.0
    RT_REC(3,2) = NYY
    RT_REC(3,3) =-NXX
    RT_REC(3,4) = U_HAT*NYY-V_HAT*NXX
    
    RT_REC(4,1) = 1.0
    RT_REC(4,2) = U_HAT+A_HAT*NXX
    RT_REC(4,3) = V_HAT+A_HAT*NYY
    RT_REC(4,4) = H_HAT+A_HAT*UC_HAT
    
    !Part 29:
    JCB_REC(:,:) = 0
    
    !Part 30:
    DO IREC=1,4
        DO JREC=1,4
            DO KREC=1,4
                
                JCB_REC(IREC,JREC) = JCB_REC(IREC,JREC)+R_REC(IREC,KREC)*EV_REC(KREC)*S_REC(KREC)*RT_REC(KREC,JREC)         !!!!!!
                
            END DO
        END DO
    END DO
            
   !Part 31: 
    DSP_REC(:) = 0
    
    !Part 32:
    DO IREC=1,4
        DO JREC=1,4
            
            DSP_REC(IREC) = DSP_REC(IREC)+JCB_REC(IREC,JREC)*DV_REC(JREC)
            
        END DO
    END DO
   
   !Part 33: 
    F1 = F1_CNV-0.5*DSP_REC(1)
    F2 = F2_CNV-0.5*DSP_REC(2)
    F3 = F3_CNV-0.5*DSP_REC(3)
    F4 = F4_CNV-0.5*DSP_REC(4)
    
   !Part 34: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
    
   !Part 35: 
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
End SUBROUTINE
    
    
!Part 36:    
Real(8) FUNCTION LOG_AVE(PRM_L,PRM_R)
   
   IMPLICIT NONE
   INTENT(IN   )::PRM_L,PRM_R
   Real(8)::PRM_L,PRM_R
   Real(8)::EPS,ZETA,F0,U0,F1
   
   EPS = 0.01
   
   ZETA = PRM_L/PRM_R
   
   F0 = (ZETA-1.0)/(ZETA+1.0)
   
   U0 = F0**2
   
   IF (U0<EPS) THEN
       
       F1 = 1+(U0/3)+((U0**2)/5)+((U0**3)/7)
       
   ELSE
       
       F1 = LOG(ZETA)/(2.0*F0)
       
   END IF
   
   LOG_AVE = (PRM_L+PRM_R)/(2.0*F1)
   
    
   RETURN
   
END FUNCTION LOG_AVE
!###########################################################################################
