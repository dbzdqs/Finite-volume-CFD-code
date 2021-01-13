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
!// Developed by: S. Kavoosi, Mechanical Eng., Amirkabir University of Technology          //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 SUBROUTINE ConMeanFlow_MatrixDiss(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,K2,K4,WNP1,WB,P,Con)
 IMPLICIT NONE
!*********************************************************************************************
 INTENT(IN   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,K2,K4,WNP1,WB,P
 INTENT(OUT  )::Con

 INTEGER::DIM,I,NC,NF1,NF2,NF,ME,NE,L,R
 Real(8)::U,V,F1,F2,F3,F4,Q,RO,RU,RV,RH,GM,A_L,A_R,M_L,M_R,M_PLUS,P_PLUS,M_MINUS,P_MINUS,&
          MM,PM,NXX,NYY,DAA
 Real(8),DIMENSION(1:4,1:DIM)::WNP1,CON
 Real(8),DIMENSION(1:5,1:DIM)::WB
 Real(8),DIMENSION(1:DIM)::P,NX,NY,DA
 INTEGER,DIMENSION(1:4,1:DIM)::IDS
 Real(8)::R_L,U_L,V_L,E_L,P_L,R_R,U_R,V_R,E_R,P_R      
 Real(8)::R_ME,U_ME,V_ME,E_ME,P_ME,RM,E
 INTEGER::II
 Real(8)::H_R,H_L
 Real(8)::UC_L,UC_R
 Real(8)::F1_L,F2_L,F3_L,F4_L
 Real(8)::F1_R,F2_R,F3_R,F4_R
 Real(8)::R_BAR,U_BAR,V_BAR,P_BAR,A_BAR,H_BAR,E_BAR
 Real(8)::K2,K4,EPS2,EPS4
 Real(8),DIMENSION(1:5,1:DIM)::LPLC,D_JST
 Real(8),DIMENSION(1:DIM):: EGV,P_PLS
 Real(8),DIMENSION(1:DIM)::P_MNS,G_EPS,G_EP1
 Real(8)::UC_BAR,JST_COEF1,JST_COEF2,JST_COEF3
 Real(8),DIMENSION(1:4,1:4)::REV_JST,LEV_JST,JCB_JST
 Real(8),DIMENSION(1:4)::EV_JST,DEL_W,ADW_JST,DEL_LW,ADLW_JST
 INTEGER::IJST,JJST,KJST
 Real(8)::KN_JST,KL_JST,R_JCB,AM,A_ME
!*********************************************************************************************	
!PART 1:
 DO I=1,NC
    CON(1,I) = 0.0
    CON(2,I) = 0.0
    CON(3,I) = 0.0
    CON(4,I) = 0.0
 END DO

!PART 2:
 DO I=NF2+1,NF

   !PART 3:
    ME = IDS(1,I)

   !PART 4:
    RM = WB(1,I)
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    E = WB(4,I)/WB(1,I)
    
   !PART 5:
    PM = WB(5,I)
    Q  = U*NX(I) + V*NY(I)
    
   !PART 6:
    F1 = Q * WB(1,I) 
    F2 = Q * WB(2,I) + PM*NX(I)
    F3 = Q * WB(3,I) + PM*NY(I)
    F4 = Q *(WB(4,I)+PM)

   !PART 7:
    CON(1,ME) = CON(1,ME) + F1
    CON(2,ME) = CON(2,ME) + F2
    CON(3,ME) = CON(3,ME) + F3
    CON(4,ME) = CON(4,ME) + F4
    
  END DO


!PART 8:
 KN_JST = 0.25                 
 KL_JST = 0.25
 
!PART 9:
 DO I=1,NC
     
    LPLC(1,I) = 0.0
    LPLC(2,I) = 0.0
    LPLC(3,I) = 0.0
    LPLC(4,I) = 0.0
    
    P_PLS(I) = 0.0
    P_MNS(I) = 0.0
    
    G_EPS(I) = 0.0
    
    D_JST(1,I) = 0.0
    D_JST(2,I) = 0.0
    D_JST(3,I) = 0.0
    D_JST(4,I) = 0.0
    
 END DO

!PART 10:
 DO I=NF1+1,NF2
     
     !PART 11:
     L = IDS(1,I)
     R = IDS(2,I)
     
     !PART 12:
     P_L = P(L)
     P_R = P(R)
     
     !PART 13:
     LPLC(1,L) = LPLC(1,L) + WNP1(1,R)-WNP1(1,L)
     LPLC(2,L) = LPLC(2,L) + WNP1(2,R)-WNP1(2,L)
     LPLC(3,L) = LPLC(3,L) + WNP1(3,R)-WNP1(3,L)
     LPLC(4,L) = LPLC(4,L) + WNP1(4,R)-WNP1(4,L)
     
     !PART 14:
     LPLC(1,R) = LPLC(1,R) + WNP1(1,L)-WNP1(1,R)
     LPLC(2,R) = LPLC(2,R) + WNP1(2,L)-WNP1(2,R)
     LPLC(3,R) = LPLC(3,R) + WNP1(3,L)-WNP1(3,R)
     LPLC(4,R) = LPLC(4,R) + WNP1(4,L)-WNP1(4,R)
     
     !PART 15:
     P_PLS(L) = P_PLS(L) + (P_R+P_L)
     P_MNS(L) = P_MNS(L) + (P_R-P_L)
     
     !PART 16:
     P_PLS(R) = P_PLS(R) + (P_L+P_R)
     P_MNS(R) = P_MNS(R) + (P_L-P_R)
     
 END DO
 
!PART 17: 
 DO I=1,NC
     
     !PART 18:    
     G_EPS(I) = DABS(P_MNS(I))/P_PLS(I)
     
 END DO
 
!PART 19:
DO I=NF1+1,NF2
     
    !PART 20:
     L = IDS(1,I)
     R = IDS(2,I)
     
     !PART 21:
     DAA = DA(I)
     NXX = NX(I)/DAA
     NYY = NY(I)/DAA
     
     !PART 22:
     R_L = WNP1(1,L)
     U_L = WNP1(2,L)/WNP1(1,L)
     V_L = WNP1(3,L)/WNP1(1,L)
     E_L = WNP1(4,L)/WNP1(1,L)
     P_L = P(L)
     
     R_R = WNP1(1,R)
     U_R = WNP1(2,R)/WNP1(1,R)
     V_R = WNP1(3,R)/WNP1(1,R)
     E_R = WNP1(4,R)/WNP1(1,R)
     P_R = P(R)
     
     !PART 23:
     A_L  = DSQRT(GM*P_L/R_L)
     A_R  = DSQRT(GM*P_R/R_R)
     
     H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
     H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
     
     !PART 24:
     EPS2 = K2*DMAX1(G_EPS(L),G_EPS(R))
     EPS4 = DMAX1( 0.D0 , (K4-EPS2) )
     
     !PART 25:
     R_BAR = 0.5*(R_L+R_R)
     U_BAR = 0.5*(U_L+U_R)
     V_BAR = 0.5*(V_L+V_R)
     E_BAR = 0.5*(E_L+E_R)
     P_BAR = 0.5*(P_L+P_R)
     A_BAR = 0.5*(A_L+A_R)
     H_BAR = 0.5*(H_L+H_R)
     
     !PART 26:
     UC_BAR = U_BAR*NXX+V_BAR*NYY

     !PART 27:
     JST_COEF1 = R_BAR/(DSQRT(2.D0)*A_BAR)
     JST_COEF2 = 1/(DSQRT(2.D0)*R_BAR*A_BAR)
     JST_COEF3 = 0.5*(GM-1)*(U_BAR*U_BAR+V_BAR*V_BAR)
       
     !PART 28:
     REV_JST(1,1) = 1
     REV_JST(1,2) = 0
     REV_JST(1,3) = JST_COEF1
     REV_JST(1,4) = JST_COEF1
     
     REV_JST(2,1) = U_BAR
     REV_JST(2,2) = R_BAR*NYY
     REV_JST(2,3) = JST_COEF1*(U_BAR+A_BAR*NXX)
     REV_JST(2,4) = JST_COEF1*(U_BAR-A_BAR*NXX)
     
     REV_JST(3,1) = V_BAR
     REV_JST(3,2) = -R_BAR*NXX
     REV_JST(3,3) = JST_COEF1*(V_BAR+A_BAR*NYY)
     REV_JST(3,4) = JST_COEF1*(V_BAR-A_BAR*NYY)
     
     REV_JST(4,1) = (JST_COEF3)/(GM-1)
     REV_JST(4,2) = R_BAR*(U_BAR*NYY-V_BAR*NXX)
     REV_JST(4,3) = JST_COEF1*((JST_COEF3+A_BAR*A_BAR)/(GM-1)+A_BAR*UC_BAR)
     REV_JST(4,4) = JST_COEF1*((JST_COEF3+A_BAR*A_BAR)/(GM-1)-A_BAR*UC_BAR)
     
     !PART 29:    
     LEV_JST(1,1) = 1-(JST_COEF3/(A_BAR*A_BAR))
     LEV_JST(1,2) = ((GM-1)*U_BAR)/(A_BAR*A_BAR)
     LEV_JST(1,3) = ((GM-1)*V_BAR)/(A_BAR*A_BAR)
     LEV_JST(1,4) = -(GM-1)/(A_BAR*A_BAR)
     
     LEV_JST(2,1) = -(U_BAR*NYY-V_BAR*NXX)/R_BAR
     LEV_JST(2,2) = NYY/R_BAR
     LEV_JST(2,3) = -NXX/R_BAR
     LEV_JST(2,4) = 0
     
     LEV_JST(3,1) = JST_COEF2*(JST_COEF3-A_BAR*UC_BAR)
     LEV_JST(3,2) = JST_COEF2*(A_BAR*NXX-(GM-1)*U_BAR)
     LEV_JST(3,3) = JST_COEF2*(A_BAR*NYY-(GM-1)*V_BAR)
     LEV_JST(3,4) = JST_COEF2*(GM-1)
     
     LEV_JST(4,1) = JST_COEF2*(JST_COEF3+A_BAR*UC_BAR)
     LEV_JST(4,2) = -JST_COEF2*(A_BAR*NXX+(GM-1)*U_BAR)
     LEV_JST(4,3) = -JST_COEF2*(A_BAR*NYY+(GM-1)*V_BAR)
     LEV_JST(4,4) = JST_COEF2*(GM-1)
     
     !PART 30:
     EV_JST(1) = DABS(UC_BAR)
     EV_JST(2) = DABS(UC_BAR)
     EV_JST(3) = DABS(UC_BAR+A_BAR*DSQRT(NXX*NXX+NYY*NYY))  
     EV_JST(4) = DABS(UC_BAR-A_BAR*DSQRT(NXX*NXX+NYY*NYY))  
     
     !PART 31:
     R_JCB = DABS(UC_BAR)+A_BAR*DSQRT(NXX*NXX+NYY*NYY)
     
     !PART 32:
     EV_JST(1) = DMAX1(EV_JST(1) , KN_JST*R_JCB)
     EV_JST(2) = DMAX1(EV_JST(2) , KN_JST*R_JCB)
     EV_JST(3) = DMAX1(EV_JST(3) , KL_JST*R_JCB)
     EV_JST(4) = DMAX1(EV_JST(4) , KL_JST*R_JCB)
     
     !PART 33:
     JCB_JST(:,:) = 0
     
     !PART 34:
     DO IJST=1,4
         DO JJST=1,4
             DO KJST=1,4
                 
                 JCB_JST(IJST,JJST) = JCB_JST(IJST,JJST)+REV_JST(IJST,KJST)*EV_JST(KJST)*LEV_JST(KJST,JJST)   
                 
             END DO
         END DO
     END DO
     
     !PART 35:
     DEL_W(1) = EPS2*(WNP1(1,R)-WNP1(1,L)) - EPS4*(LPLC(1,R)-LPLC(1,L))
     DEL_W(2) = EPS2*(WNP1(2,R)-WNP1(2,L)) - EPS4*(LPLC(2,R)-LPLC(2,L))
     DEL_W(3) = EPS2*(WNP1(3,R)-WNP1(3,L)) - EPS4*(LPLC(3,R)-LPLC(3,L))
     DEL_W(4) = EPS2*(WNP1(4,R)-WNP1(4,L)) - EPS4*(LPLC(4,R)-LPLC(4,L))
     
     !PART 36: 
     ADW_JST(:) = 0
     
     !PART 37:
     DO IJST=1,4
         DO JJST=1,4
             
             ADW_JST(IJST) = ADW_JST(IJST)+JCB_JST(IJST,JJST)*DEL_W(JJST)
             
         END DO
     END DO
     
     !PART 38:
     D_JST(1,L) = D_JST(1,L) + ADW_JST(1)*DAA
     D_JST(2,L) = D_JST(2,L) + ADW_JST(2)*DAA
     D_JST(3,L) = D_JST(3,L) + ADW_JST(3)*DAA
     D_JST(4,L) = D_JST(4,L) + ADW_JST(4)*DAA
     
     !PART 39:
     D_JST(1,R) = D_JST(1,R) - ADW_JST(1)*DAA
     D_JST(2,R) = D_JST(2,R) - ADW_JST(2)*DAA
     D_JST(3,R) = D_JST(3,R) - ADW_JST(3)*DAA
     D_JST(4,R) = D_JST(4,R) - ADW_JST(4)*DAA
     
END DO

!PART 40:
 DO I=NF1+1,NF2

   !PART 41:
   L = IDS(1,I)
   R = IDS(2,I)

   !PART 42:
   DAA = DA(I)
   NXX = NX(I)/DAA
   NYY = NY(I)/DAA

   !PART 43:
   R_L = WNP1(1,L)
   U_L = WNP1(2,L)/WNP1(1,L)
   V_L = WNP1(3,L)/WNP1(1,L)
   E_L = WNP1(4,L)/WNP1(1,L)
   P_L = P(L)
   
   R_R = WNP1(1,R)
   U_R = WNP1(2,R)/WNP1(1,R)
   V_R = WNP1(3,R)/WNP1(1,R)
   E_R = WNP1(4,R)/WNP1(1,R)
   P_R = P(R)
   
   !PART 44: 
   H_L = (WNP1(4,L)+P(L))/WNP1(1,L)
   H_R = (WNP1(4,R)+P(R))/WNP1(1,R)
   
   A_L  = DSQRT(GM*P(L)/WNP1(1,L))
   A_R  = DSQRT(GM*P(R)/WNP1(1,R))
   
   !PART 45:
   UC_L = U_L*NXX+V_L*NYY
   UC_R = U_R*NXX+V_R*NYY
   
   !PART 46:
   F1_L = UC_L*R_L
   F2_L = UC_L*WNP1(2,L)+P(L)*NXX
   F3_L = UC_L*WNP1(3,L)+P(L)*NYY
   F4_L = UC_L*(WNP1(4,L)+P(L))
   
   F1_R = UC_R*R_R
   F2_R = UC_R*WNP1(2,R)+P(R)*NXX
   F3_R = UC_R*WNP1(3,R)+P(R)*NYY
   F4_R = UC_R*(WNP1(4,R)+P(R))
   
   !PART 47:
   F1 = 0.5*(F1_L+F1_R)
   F2 = 0.5*(F2_L+F2_R)
   F3 = 0.5*(F3_L+F3_R)
   F4 = 0.5*(F4_L+F4_R)
      
   !PART 48: 
   Con(1,L) = Con(1,L) + F1*DAA 
   Con(2,L) = Con(2,L) + F2*DAA 
   Con(3,L) = Con(3,L) + F3*DAA 
   Con(4,L) = Con(4,L) + F4*DAA 
                                
   !PART 49:                     
   Con(1,R) = Con(1,R) - F1*DAA 
   Con(2,R) = Con(2,R) - F2*DAA 
   Con(3,R) = Con(3,R) - F3*DAA 
   Con(4,R) = Con(4,R) - F4*DAA 
   
 End Do

 !PART 50:
 Con(1,:) = Con(1,:) - D_JST(1,:)
 Con(2,:) = Con(2,:) - D_JST(2,:)
 Con(3,:) = Con(3,:) - D_JST(3,:)
 Con(4,:) = Con(4,:) - D_JST(4,:)
!*********************************************************************************************
 End
!###########################################################################################