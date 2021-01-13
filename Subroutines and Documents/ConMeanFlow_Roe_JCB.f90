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
 Subroutine ConMeanFlow_Roe_JCB(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R,IROE,JROE,KROE,IENT
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA,R_L,U_L,V_L,H_L,R_R,U_R,V_R,H_R,R_RT,R_ROE,U_ROE,V_ROE,A_ROE,&
		  H_ROE,UC_ROE,UC_L,UC_R,F1_L,F2_L,F3_L,F4_L,F1_R,F2_R,F3_R,F4_R,ROE_COEF1,ROE_COEF2,&
		  ROE_COEF3,EV1_R,EV2_R,EV3_R,EV1_L,EV2_L,EV3_L,EV1_ROE,EV2_ROE,EV3_ROE,EPS1,EPS2,&
		  EPS3,EV1_NEW,EV2_NEW,EV3_NEW
 Real(8),Dimension(1:4,1:Dim)::WNP1,Con
 Real(8),Dimension(1:5,1:Dim)::WB
 Real(8),Dimension(1:Dim)::P,NX,NY,DA
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),DIMENSION(1:4)::DEL_W,ADW_ROE,EV_ROE
 Real(8),DIMENSION(1:4,1:4)::REV_ROE,LEV_ROE,JCB_ROE
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
    ROE_COEF1 = R_ROE/(DSQRT(2.D0)*A_ROE)
    ROE_COEF2 = 1/(DSQRT(2.D0)*R_ROE*A_ROE)
    ROE_COEF3 = 0.5*(GM-1)*(U_ROE*U_ROE+V_ROE*V_ROE)
    
   !Part 17:    
    REV_ROE(1,1) = 1
    REV_ROE(1,2) = 0
    REV_ROE(1,3) = ROE_COEF1
    REV_ROE(1,4) = ROE_COEF1
    
    REV_ROE(2,1) = U_ROE
    REV_ROE(2,2) = R_ROE*NYY
    REV_ROE(2,3) = ROE_COEF1*(U_ROE+A_ROE*NXX)
    REV_ROE(2,4) = ROE_COEF1*(U_ROE-A_ROE*NXX)
    
    REV_ROE(3,1) = V_ROE
    REV_ROE(3,2) = -R_ROE*NXX
    REV_ROE(3,3) = ROE_COEF1*(V_ROE+A_ROE*NYY)
    REV_ROE(3,4) = ROE_COEF1*(V_ROE-A_ROE*NYY)
    
    REV_ROE(4,1) = (ROE_COEF3)/(GM-1)
    REV_ROE(4,2) = R_ROE*(U_ROE*NYY-V_ROE*NXX)
    REV_ROE(4,3) = ROE_COEF1*((ROE_COEF3+A_ROE*A_ROE)/(GM-1)+A_ROE*UC_ROE)
    REV_ROE(4,4) = ROE_COEF1*((ROE_COEF3+A_ROE*A_ROE)/(GM-1)-A_ROE*UC_ROE)
    
    
   !Part 18:    
    LEV_ROE(1,1) = 1-(ROE_COEF3/(A_ROE*A_ROE))
    LEV_ROE(1,2) = ((GM-1)*U_ROE)/(A_ROE*A_ROE)
    LEV_ROE(1,3) = ((GM-1)*V_ROE)/(A_ROE*A_ROE)
    LEV_ROE(1,4) = -(GM-1)/(A_ROE*A_ROE)
    
    LEV_ROE(2,1) = -(U_ROE*NYY-V_ROE*NXX)/R_ROE
    LEV_ROE(2,2) = NYY/R_ROE
    LEV_ROE(2,3) = -NXX/R_ROE
    LEV_ROE(2,4) = 0
    
    LEV_ROE(3,1) = ROE_COEF2*(ROE_COEF3-A_ROE*UC_ROE)
    LEV_ROE(3,2) = ROE_COEF2*(A_ROE*NXX-(GM-1)*U_ROE)
    LEV_ROE(3,3) = ROE_COEF2*(A_ROE*NYY-(GM-1)*V_ROE)
    LEV_ROE(3,4) = ROE_COEF2*(GM-1)
    
    LEV_ROE(4,1) = ROE_COEF2*(ROE_COEF3+A_ROE*UC_ROE)
    LEV_ROE(4,2) = -ROE_COEF2*(A_ROE*NXX+(GM-1)*U_ROE)
    LEV_ROE(4,3) = -ROE_COEF2*(A_ROE*NYY+(GM-1)*V_ROE)
    LEV_ROE(4,4) = ROE_COEF2*(GM-1)
    
    
    
   !Part 19:
    IENT = 2
    !!ENTROPY FIX CORRECTIONS!
    IF (IENT /= 0) THEN
        
        !Part 20:
        EV1_R = UC_R
        EV2_R = UC_R+A_R*DSQRT(NXX*NXX+NYY*NYY)
        EV3_R = UC_R-A_R*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_L = UC_L
        EV2_L = UC_L+A_L*DSQRT(NXX*NXX+NYY*NYY)
        EV3_L = UC_L-A_L*DSQRT(NXX*NXX+NYY*NYY)
        
        EV1_ROE = UC_ROE
        EV2_ROE = UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        EV3_ROE = UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY)
        
        !Part 21:
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
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 22:
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
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 23:    
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
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 24:
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
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        !Part 25:
        ELSE IF (IENT == 5) THEN
            
            !!ENTROPY FIX HOFFMAN-CHIANG!
            EPS1 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            EPS2 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            EPS3 = 0.1           !EPS SHOULD BE: 0<EPS<0.125
            
            
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
            
            
            EV_ROE(1) = DABS(EV1_NEW)
            EV_ROE(2) = DABS(EV1_NEW)
            EV_ROE(3) = DABS(EV2_NEW)  
            EV_ROE(4) = DABS(EV3_NEW)
        
        
        END IF
        
        
    ELSE
        !Part 26:
        EV_ROE(1) = DABS(UC_ROE)
        EV_ROE(2) = DABS(UC_ROE)
        EV_ROE(3) = DABS(UC_ROE+A_ROE*DSQRT(NXX*NXX+NYY*NYY))  
        EV_ROE(4) = DABS(UC_ROE-A_ROE*DSQRT(NXX*NXX+NYY*NYY))  
    
    END IF
    
   
        
    
   !Part 27:    
    JCB_ROE(:,:) = 0
    
    !Part 28:
    DO IROE=1,4
        DO JROE=1,4
            DO KROE=1,4
                
                JCB_ROE(IROE,JROE) = JCB_ROE(IROE,JROE)+REV_ROE(IROE,KROE)*EV_ROE(KROE)*LEV_ROE(KROE,JROE)         !!!!!!
                
            END DO
        END DO
    END DO
            
    
   !Part 29:
    DEL_W(1) = WNP1(1,R)-WNP1(1,L)
    DEL_W(2) = WNP1(2,R)-WNP1(2,L)
    DEL_W(3) = WNP1(3,R)-WNP1(3,L)
    DEL_W(4) = WNP1(4,R)-WNP1(4,L)
   
    
   !Part 30:    
    ADW_ROE(:) = 0
    
    !Part 31:
    DO IROE=1,4
        DO JROE=1,4
            
            ADW_ROE(IROE) = ADW_ROE(IROE)+JCB_ROE(IROE,JROE)*DEL_W(JROE)
            
        END DO
    END DO
    
   !Part 32:
    F1 = 0.5*(F1_L+F1_R-ADW_ROE(1))*DAA
    F2 = 0.5*(F2_L+F2_R-ADW_ROE(2))*DAA
    F3 = 0.5*(F3_L+F3_R-ADW_ROE(3))*DAA
    F4 = 0.5*(F4_L+F4_R-ADW_ROE(4))*DAA
    
   !Part 33:    
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4
    
   !Part 34:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
