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
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_ScalarDiss_ALE_3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,GF,Con)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,WB,P,GF
 Intent (Out  )::Con

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::U,V,W,NXX,NYY,NZZ,F1,F2,F3,F4,F5,Q,Pm,R,RU,RV,RW,RE,GM,K2,K4
 Real(8),Dimension(1:5,1:Dim)::Con,Dis,WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P,GF
 Integer,Dimension(1:6,1:Dim)::IDS
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    Con(1,I) = 0.0
    Con(2,I) = 0.0
    Con(3,I) = 0.0
    Con(4,I) = 0.0
    Con(5,I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF
 
   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)
    
    NXX = NX(I)    
	NYY = NY(I)    
	NZZ = NZ(I)

   !Part 5:
    Q  = U*NXX + V*NYY + W*NZZ
    Pm = WB(6,I)

   !Part 6:
    F1 = (Q - GF(I) ) *  Wb(1,I) 
    F2 = (Q - GF(I) ) *  Wb(2,I) + Pm*NXX
    F3 = (Q - GF(I) ) *  Wb(3,I) + Pm*NYY
    F4 = (Q - GF(I) ) *  Wb(4,I) + Pm*NZZ
    F5 = (Q - GF(I) ) *  Wb(5,I) + Pm* Q 

   !Part 7:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
    Con(5,ME) = Con(5,ME) + F5
	
 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 10:
    NXX = NX(I)
    NYY = NY(I)
    NZZ = Nz(I)

   !Part 11:
    R  = 0.5 * ( WNP1(1,ME) + WNP1(1,NE) )
	RU = 0.5 * ( WNP1(2,ME) + WNP1(2,NE) )
	RV = 0.5 * ( WNP1(3,ME) + WNP1(3,NE) )
	RW = 0.5 * ( WNP1(4,ME) + WNP1(4,NE) )
    RE = 0.5 * ( WNP1(5,ME) + WNP1(5,NE) )
    Pm = 0.5 * ( P(ME)      + P(NE)      )

   !Part 12:
    U = RU / R
    V = RV / R
    W = RW / R

   !Part 13:
    Q = U*NXX + V*NYY + W*NZZ

   !Part 14:
    F1 = (Q - GF(I) ) * R 
    F2 = (Q - GF(I) ) * RU + Pm*NXX
    F3 = (Q - GF(I) ) * RV + Pm*NYY
    F4 = (Q - GF(I) ) * RW + Pm*NZZ
    F5 = (Q - GF(I) ) * RE + Pm*Q 

   !Part 15:
    Con(1,ME) = Con(1,ME) + F1
    Con(2,ME) = Con(2,ME) + F2
    Con(3,ME) = Con(3,ME) + F3
    Con(4,ME) = Con(4,ME) + F4
    Con(5,ME) = Con(5,ME) + F5

   !Part 16:
    Con(1,NE) = Con(1,NE) - F1
    Con(2,NE) = Con(2,NE) - F2
    Con(3,NE) = Con(3,NE) - F3
    Con(4,NE) = Con(4,NE) - F4
    Con(5,NE) = Con(5,NE) - F5
 
 End Do

!Part 17:
 Call DisJameson_ALE_3D(Dim,NC,NF1,NF2,IDS,NX,NY,NZ,DA,GM,K2,K4,WNP1,P,GF,Dis)

!Part 18:
 DO I=1,NC
    Con(1,I) = Con(1,I) - Dis(1,I)
    Con(2,I) = Con(2,I) - Dis(2,I)
    Con(3,I) = Con(3,I) - Dis(3,I)
    Con(4,I) = Con(4,I) - Dis(4,I)
    Con(5,I) = Con(5,I) - Dis(5,I)
 End Do
!*********************************************************************************************
 End
!###########################################################################################

