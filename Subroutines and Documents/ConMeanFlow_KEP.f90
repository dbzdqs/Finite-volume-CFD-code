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
 Subroutine ConMeanFlow_KEP(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
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
 Real(8):: R_BAR,A_BAR,P_BAR,U_BAR,V_BAR,E_BAR,H_BAR,UC_BAR
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
   
   !Part 13: 
    R_BAR = 0.5*(R_L+R_R)
    U_BAR = 0.5*(U_L+U_R)
    V_BAR = 0.5*(V_L+V_R)
    E_BAR = 0.5*(E_L+E_R)
    P_BAR = 0.5*(P_L+P_R)
    H_BAR = 0.5*(H_L+H_R)
    
   !Part 14: 
    UC_BAR = U_BAR*NXX+V_BAR*NYY
    
   !Part 15:
    F1 = (UC_BAR*R_BAR               )
    F2 = (UC_BAR*R_BAR*U_BAR+P_BAR*NXX)
    F3 = (UC_BAR*R_BAR*V_BAR+P_BAR*NYY)
    F4 = (UC_BAR*R_BAR*H_BAR )
    
   !Part 16: 
    Con(1,L) = Con(1,L) + F1*DAA
    Con(2,L) = Con(2,L) + F2*DAA
    Con(3,L) = Con(3,L) + F3*DAA
    Con(4,L) = Con(4,L) + F4*DAA
                            
   !Part 17:                
    Con(1,R) = Con(1,R) - F1*DAA
    Con(2,R) = Con(2,R) - F2*DAA
    Con(3,R) = Con(3,R) - F3*DAA
    Con(4,R) = Con(4,R) - F4*DAA
    
 End Do
!*********************************************************************************************
 End SUBROUTINE
!###########################################################################################
