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
!// Developed by: S. Kavoosi, Mechanical Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_Roe_ENT(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R,IENT
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA,R_L,U_L,V_L,H_L,R_R,U_R,V_R,H_R,R_RT,R_ROE,U_ROE,V_ROE,A_ROE,&
		  H_ROE,UC_ROE,DEL_UC,UC_L,UC_R,F1_L,F2_L,F3_L,F4_L,F1_R,F2_R,F3_R,F4_R,DEL_R,DEL_U,&
		  DEL_V,DEL_P,ROE_COEF1,ROE_COEF2,ROE_COEF3,ROE_COEF4,ROE_COEF5,F1_ROE1,F2_ROE1,&
		  F3_ROE1,F4_ROE1,F1_ROE4,F2_ROE4,F3_ROE4,F4_ROE4,F1_ROE5,F2_ROE5,F3_ROE5,F4_ROE5,&
		  EV1_R,EV2_R,EV3_R,EV1_L,EV2_L,EV3_L,EV1_ROE,EV2_ROE,EV3_ROE,EPS1,EPS2,EPS3,EV1_NEW,&
		  EV2_NEW,EV3_NEW
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
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
    
    R_R = WNP1(1,R)
    U_R = WNP1(2,R)/WNP1(1,R)
    V_R = WNP1(3,R)/WNP1(1,R)
    
   !Part 12: 
    H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
    H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
    
    A_L  = DSQRT(GM*P(L)/WNP1(1,L))
	A_R  = DSQRT(GM*P(R)/WNP1(1,R))
    
   !Part 13: 
    R_RT = DSQRT(R_R/R_L)

    R_ROE = DSQRT(R_L*R_R)
    U_ROE = (U_L+U_R*R_RT)/(1.0+R_RT)
    V_ROE = (V_L+V_R*R_RT)/(1.0+R_RT)
    H_ROE = (H_L+H_R*R_RT)/(1.0+R_RT)
    
    A_ROE = DSQRT((GM-1.0)*(H_ROE-0.5*(U_ROE*U_ROE+V_ROE*V_ROE)))
    
   !Part 14: 
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
    UC_ROE = U_ROE*NXX+V_ROE*NYY
    
    DEL_UC = (U_R-U_L)*NXX+(V_R-V_L)*NYY
    
   !Part 15: 
    F1_L = UC_L*R_L
    F2_L = UC_L*WNP1(2,L)+P(L)*NXX
    F3_L = UC_L*WNP1(3,L)+P(L)*NYY
    F4_L = UC_L*(WNP1(4,L)+P(L))
    
    F1_R = UC_R*R_R
    F2_R = UC_R*WNP1(2,R)+P(R)*NXX
    F3_R = UC_R*WNP1(3,R)+P(R)*NYY
    F4_R = UC_R*(WNP1(4,R)+P(R))
    
   !Part 16: 
    DEL_R = R_R-R_L
    DEL_U = U_R-U_L
    DEL_V = V_R-V_L
    DEL_P = P(R)-P(L)
    
    
   !Part 17: 
    IENT = 2
    !!ENTROPY FIX CORRECTIONS!
    IF (IENT /= 0) THEN
        
        !Part 18:
        EV1_R = UC_R
        EV2_R = UC_R+A_R*DSQRT(NXX*NXX+NYY*NYY)
        EV3_R = UC_R-A_R*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_L = UC_L
        EV2_L = UC_L+A_L*DSQRT(NXX*NXX+NYY*NYY)
        EV3_L = UC_L-A_L*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_ROE = UC_ROE
        EV2_ROE = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        EV3_ROE = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        
        !Part 19:
        IF (IENT == 1) THEN
            
            !!ENTROPY FIX KERMANI#1!
            EPS1 = 2.D0*DMAX1(0.D0,(EV1_R-EV1_L))
            EPS2 = 2.D0*DMAX1(0.D0,(EV2_R-EV2_L))
            EPS3 = 2.D0*DMAX1(0.D0,(EV3_R-EV3_L))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 20:    
        ELSE IF (IENT == 2) THEN
            
            !!ENTROPY FIX KERMANI#2!
            EPS1 = 4.D0*DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = 4.D0*DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = 4.D0*DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 21:    
        ELSE IF (IENT == 3) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#1!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = (EV1_ROE*EV1_ROE+EPS1*EPS1)/(2*EPS1)
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = (EV2_ROE*EV2_ROE+EPS2*EPS2)/(2*EPS2)
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = (EV3_ROE*EV3_ROE+EPS3*EPS3)/(2*EPS3)
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 22:    
        ELSE IF (IENT == 4) THEN
            
            !!ENTROPY FIX HARTEN-HYMAN#2!
            EPS1 = DMAX1(0.D0,(EV1_ROE-EV1_L),(EV1_R-EV1_ROE))
            EPS2 = DMAX1(0.D0,(EV2_ROE-EV2_L),(EV2_R-EV2_ROE))
            EPS3 = DMAX1(0.D0,(EV3_ROE-EV3_L),(EV3_R-EV3_ROE))
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
        !Part 23:
        ELSE IF (IENT == 5) THEN
            
            !!ENTROPY FIX HOFFMAN-CHIANG!
            EPS1 = 0.05          !EPS SHOULD BE: 0<EPS<0.125
            EPS2 = 0.05           !EPS SHOULD BE: 0<EPS<0.125
            EPS3 = 0.05           !EPS SHOULD BE: 0<EPS<0.125
            
            
            IF (DABS(EV1_ROE)<EPS1) THEN
                
                EV1_NEW = EPS1
                
            ELSE
                
                EV1_NEW = EV1_ROE
                
            END IF
            
            IF (DABS(EV2_ROE)<EPS2) THEN
                
                EV2_NEW = EPS2
                
            ELSE
                
                EV2_NEW = EV2_ROE
                
            END IF
            
            IF (DABS(EV3_ROE)<EPS3) THEN
                
                EV3_NEW = EPS3
                
            ELSE
                
                EV3_NEW = EV3_ROE
                
            END IF
            
            
        END IF
        
        
    ELSE
        !Part 24:
        EV1_NEW = UC_ROE
        EV2_NEW = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)  
        EV3_NEW = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)  
    
    END IF
    
    
   !Part 25: 
    ROE_COEF1 = DABS(EV1_NEW)*(DEL_R-DEL_P/(A_ROE*A_ROE))
    ROE_COEF2 = DABS(EV1_NEW)*R_ROE
    ROE_COEF3 = 0.5*(U_ROE*U_ROE+V_ROE*V_ROE)
    ROE_COEF4 = DABS(EV2_NEW)*((DEL_P+(R_ROE*A_ROE*DEL_UC))/(2.0*A_ROE*A_ROE))
    ROE_COEF5 = DABS(EV3_NEW)*((DEL_P-(R_ROE*A_ROE*DEL_UC))/(2.0*A_ROE*A_ROE))
    
   !Part 26: 
    F1_ROE1 = ROE_COEF1*1.0+ROE_COEF2*0.0
    F2_ROE1 = ROE_COEF1*U_ROE+ROE_COEF2*(DEL_U-NXX*DEL_UC)
    F3_ROE1 = ROE_COEF1*V_ROE+ROE_COEF2*(DEL_V-NYY*DEL_UC)
    F4_ROE1 = ROE_COEF1*ROE_COEF3+ROE_COEF2*(U_ROE*DEL_U+V_ROE*DEL_V-UC_ROE*DEL_UC)
    
    F1_ROE4 = ROE_COEF4*(1.0)
    F2_ROE4 = ROE_COEF4*(U_ROE+A_ROE*NXX)
    F3_ROE4 = ROE_COEF4*(V_ROE+A_ROE*NYY)
    F4_ROE4 = ROE_COEF4*(H_ROE+A_ROE*UC_ROE)
    
    F1_ROE5 = ROE_COEF5*(1.0)
    F2_ROE5 = ROE_COEF5*(U_ROE-A_ROE*NXX)
    F3_ROE5 = ROE_COEF5*(V_ROE-A_ROE*NYY)
    F4_ROE5 = ROE_COEF5*(H_ROE-A_ROE*UC_ROE)
    
    
   !Part 27: 
    F1 = 0.5*(F1_L+F1_R-(F1_ROE1+F1_ROE4+F1_ROE5))
    F2 = 0.5*(F2_L+F2_R-(F2_ROE1+F2_ROE4+F2_ROE5))
    F3 = 0.5*(F3_L+F3_R-(F3_ROE1+F3_ROE4+F3_ROE5))
    F4 = 0.5*(F4_L+F4_R-(F4_ROE1+F4_ROE4+F4_ROE5))
    
   !Part 28: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
    
   !Part 29: 
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End 
!###########################################################################################
