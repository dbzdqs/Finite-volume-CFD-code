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
 Subroutine DisJameson_ALE_3D(Dim,NC,NF1,NF2,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,P,GF,Dis)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,P,GF
 Intent(Out  )::Dis

 Integer::Dim,I,NC,NF1,NF2,ME,NE
 Real(8)::C2,GM,U,V,W,Ai,E2,E4,T1,T2,T3,T4,T5,S1,S2,S3,S4,S5,F1,F2,F3,F4,F5
 Real(8)::DW1,DW2,DW3,DW4,DW5,K2,K4
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::DW,Dis,WNP1
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P,GF
!*********************************************************************************************
!Part 1:
 DO I=1,NC
    DW(1,I) = 0.0
    DW(2,I) = 0.0
    DW(3,I) = 0.0
    DW(4,I) = 0.0
    DW(5,I) = 0.0
    
    Dis(1,I) = 0.0
    Dis(2,I) = 0.0
    Dis(3,I) = 0.0
    Dis(4,I) = 0.0
    Dis(5,I) = 0.0
 End Do
 
!Part 2:
 DO I=NF1+1,NF2
     
   !Part 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
    
   !Part 4:
    DW1 = WNP1(1,NE)-WNP1(1,ME)
    DW2 = WNP1(2,NE)-WNP1(2,ME)
    DW3 = WNP1(3,NE)-WNP1(3,ME)
    DW4 = WNP1(4,NE)-WNP1(4,ME)
    DW5 = WNP1(5,NE)-WNP1(5,ME)  + P(NE)-P(ME)
    
   !Part 5:
    DW(1,ME) = DW(1,ME) + DW1
    DW(2,ME) = DW(2,ME) + DW2
    DW(3,ME) = DW(3,ME) + DW3
    DW(4,ME) = DW(4,ME) + DW4
    DW(5,ME) = DW(5,ME) + DW5
    
   !Part 6:
    DW(1,NE) = DW(1,NE) - DW1
    DW(2,NE) = DW(2,NE) - DW2
    DW(3,NE) = DW(3,NE) - DW3
    DW(4,NE) = DW(4,NE) - DW4
    DW(5,NE) = DW(5,NE) - DW5
 End Do
 
!Part 7:
 DO I=NF1+1,NF2
     
   !Part 8:
    ME = IDS(1,I)
    NE = IDS(2,I)
    
   !Part 9:
    C2 = GM*(P(ME)+P(NE))/(WNP1(1,ME)+WNP1(1,NE))
    
   !Part 10:
    U = (WNP1(2,ME)/WNP1(1,ME)+WNP1(2,NE)/WNP1(1,NE))*0.5
    V = (WNP1(3,ME)/WNP1(1,ME)+WNP1(3,NE)/WNP1(1,NE))*0.5
    W = (WNP1(4,ME)/WNP1(1,ME)+WNP1(4,NE)/WNP1(1,NE))*0.5
    
   !Part 11:
   !Ai = ABS(U*NX(I)+V*NX(I)+W*NZ(I)-GF(I)) + SQRT(C2*DA(I)*DA(I)) !????????
    Ai = ABS(U*NX(I)+V*NY(I)+W*NZ(I)-GF(I)) + SQRT(C2*DA(I)*DA(I)) !????????
    
   !Part 12:
    E2 = K2*ABS((P(NE)-P(ME))/(P(NE)+P(ME)))
    E4 = DMAX1(0.0,(K4-E2))
    
   !Part 13:
    F1 = E4*(DW(1,ME)-DW(1,NE))
    F2 = E4*(DW(2,ME)-DW(2,NE))
    F3 = E4*(DW(3,ME)-DW(3,NE))
    F4 = E4*(DW(4,ME)-DW(4,NE))
    F5 = E4*(DW(5,ME)-DW(5,NE))
    
    DW1 = WNP1(1,ME)-WNP1(1,NE)
    DW2 = WNP1(2,ME)-WNP1(2,NE)
    DW3 = WNP1(3,ME)-WNP1(3,NE)
    DW4 = WNP1(4,ME)-WNP1(4,NE)
    DW5 = WNP1(5,ME)-WNP1(5,NE) + P(ME)-P(NE)
    
    S1 = E2*DW1
    S2 = E2*DW2
    S3 = E2*DW3
    S4 = E2*DW4
    S5 = E2*DW5
    
    T1 = Ai*(F1-S1)
    T2 = Ai*(F2-S2)
    T3 = Ai*(F3-S3)
    T4 = Ai*(F4-S4)
    T5 = Ai*(F5-S5)
    
   !Part 14:
    Dis(1,ME) = Dis(1,ME) + T1
    Dis(2,ME) = Dis(2,ME) + T2
    Dis(3,ME) = Dis(3,ME) + T3
    Dis(4,ME) = Dis(4,ME) + T4
    Dis(5,ME) = Dis(5,ME) + T5
    
   !Part 15:
    Dis(1,NE) = Dis(1,NE) - T1
    Dis(2,NE) = Dis(2,NE) - T2
    Dis(3,NE) = Dis(3,NE) - T3
    Dis(4,NE) = Dis(4,NE) - T4
    Dis(5,NE) = Dis(5,NE) - T5
    
 End Do 
!*********************************************************************************************
 End
!###########################################################################################
