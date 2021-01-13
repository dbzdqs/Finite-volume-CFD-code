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
!// Date: June, 10, 2017                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_CUSP_Zha(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P
 Intent(Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,Ro,RU,RV,RH,GM,a_L,a_R,M_L,M_R,M_Plus,P_Plus,M_Minus,P_Minus,&
          Mm,Pm,Nxx,Nyy,DAA
 Real(8),Dimension(1:4,1:110000)::WNP1,Con
 Real(8),Dimension(1:5,1:110000)::WB
 Real(8),Dimension(1:110000)::P,NX,NY,DA
 Integer,Dimension(1:4,1:110000)::IDS
 
 Real(8)::R_L,U_L,V_L,H_L,P_L,E_L
 Real(8)::R_R,U_R,V_R,H_R,P_R,E_R
 Real(8)::R_RT
 Real(8)::R_ROE,U_ROE,V_ROE,A_ROE,H_ROE
 Real(8)::UC_ROE,DEL_UC,UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8):: R_BAR,A_BAR,P_BAR,U_BAR,V_BAR,E_BAR,H_BAR,UC_BAR
 Real(8)::AC_L,AC_R,AC_BAR,M_BAR
 Real(8)::DLT_P,DLT_M,M_LP,M_LM,M_RP,M_RM
 Real(8)::BET_L,BET_R,M_HALF,PHI,M_HALFP,M_HALFM
 Real(8)::ALF_PL,ALF_MR,AC_P,AC_M,P_PL,P_MR,S_PL,S_MR,D_PL,D_MR
 Real(8)::F1_C,F2_C,F3_C,F4_C
 Real(8)::F1_P,F2_P,F3_P,F4_P
 
 Real(8)::MM_L,MM_R
 Real(8)::ALF_L,ALF_R,UP_L,UM_R,RU_HALF,ALF
 INTEGER::IZHA
 
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
    
    A_L = DSQRT(GM*P(L)/WNP1(1,L))
	A_R = DSQRT(GM*P(R)/WNP1(1,R))
    
    !Part 13: 
    AC_BAR = 0.5*(A_L+A_R)
    
    !Part 14:
    UC_L = U_L*NXX+V_L*NYY
    UC_R = U_R*NXX+V_R*NYY
    
    !Part 15:
    M_L = UC_L/AC_BAR
    M_R = UC_R/AC_BAR
    M_BAR = 0.5*(M_L+M_R)
        
    !Part 16:
	IF(DABS(M_L)>1.0)THEN
	 M_PLUS = 0.5*(M_L+DABS(M_L))
	 
	ELSE
     M_PLUS = 0.25*(M_L+1.)*(M_L+1.)
	 
    END IF

    !Part 17:
	IF(DABS(M_R)>1.0)THEN
	 M_MINUS = 0.5*(M_R-DABS(M_R))
	 
	ELSE
     M_MINUS =-0.25*(M_R-1.)*(M_R-1.)
	 
    END IF

    !Part 18:
    MM = M_PLUS+M_MINUS
    
       
   !PART 19:
	IF(MM<=-1.0)THEN
        
        F1 = UC_R*R_R
        F2 = UC_R*WNP1(2,R)+P(R)*NXX
        F3 = UC_R*WNP1(3,R)+P(R)*NYY
        F4 = UC_R*(WNP1(4,R)+P(R))
        
    !Part 20:    
    ELSE IF(MM>=1.0)THEN
        
        F1 = UC_L*R_L
        F2 = UC_L*WNP1(2,L)+P(L)*NXX
        F3 = UC_L*WNP1(3,L)+P(L)*NYY
        F4 = UC_L*(WNP1(4,L)+P(L))
    
    !Part 21:
    ELSE
        
        !Part 22:
        !!ZHA-CUSP
        ALF_L = 2.0*((P_L/R_L)/((P_L/R_L)+(P_R/R_R)))
        ALF_R = 2.0*((P_R/R_R)/((P_L/R_L)+(P_R/R_R)))
        
        !Part 23:
        UP_L = AC_BAR*( 0.5*(M_L+DABS(M_L)) + ALF_L*(0.25*((M_L+1.0)**2)-0.5*(M_L+DABS(M_L))) )
        UM_R = AC_BAR*( 0.5*(M_R-DABS(M_R)) + ALF_R*(-0.25*((M_R-1.0)**2)-0.5*(M_R-DABS(M_R))) )
        
        !Part 24:
        RU_HALF = R_L*UP_L+R_R*UM_R
        
        !Part 25:
        F1_C = 0.5*(RU_HALF*(1.0+1.0)-DABS(RU_HALF)*(1.0-1.0))
        F2_C = 0.5*(RU_HALF*(U_L+U_R)-DABS(RU_HALF)*(U_R-U_L))
        F3_C = 0.5*(RU_HALF*(V_L+V_R)-DABS(RU_HALF)*(V_R-V_L))
        
        !Part 26:
        IZHA = 1
        
        !Part 27:
        IF (IZHA==1) THEN
            
            !!ZHA-CUSP2
            ALF_L = 2.0*((H_L/R_L)/((H_L/R_L)+(H_R/R_R)))
            ALF_R = 2.0*((H_R/R_R)/((H_L/R_L)+(H_R/R_R)))
            
            UP_L = AC_BAR*( 0.5*(M_L+DABS(M_L)) + ALF_L*(0.25*((M_L+1.0)**2)-0.5*(M_L+DABS(M_L))) )
            UM_R = AC_BAR*( 0.5*(M_R-DABS(M_R)) + ALF_R*(-0.25*((M_R-1.0)**2)-0.5*(M_R-DABS(M_R))) )
            
            RU_HALF = R_L*UP_L+R_R*UM_R
            
        END IF
        
        !Part 28:
        F4_C = 0.5*(RU_HALF*(E_L+E_R)-DABS(RU_HALF)*(E_R-E_L))
        
        !Part 29:
        ALF = 3.0/16.0
        P_PL = 0.25*((M_L+1.0)**2)*(2.0-M_L)+ALF*M_L*((M_L**2-1.0)**2)
        P_MR = 0.25*((M_R-1.0)**2)*(2.0+M_R)-ALF*M_R*((M_R**2-1.0)**2)
        
        !Part 30:
        F1_P = 0.0
        F2_P = P_PL*P_L*NXX+P_MR*P_R*NXX
        F3_P = P_PL*P_L*NYY+P_MR*P_R*NYY
        F4_P = 0.5*(P_L*(UC_L+AC_BAR)+P_R*(UC_R-AC_BAR))
        
        !Part 31:
        F1 = F1_C + F1_P 
        F2 = F2_C + F2_P 
        F3 = F3_C + F3_P 
        F4 = F4_C + F4_P 
                
	ENDIF

    
   !Part 32: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 33:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
