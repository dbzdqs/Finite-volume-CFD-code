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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
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
    a_L  = Dsqrt(GM*P(L)/WNP1(1,L))
	a_R  = Dsqrt(GM*P(R)/WNP1(1,R))

	M_L = (WNP1(2,L)*NXX+WNP1(3,L)*NYY) / (WNP1(1,L) * a_L)
	M_R = (WNP1(2,R)*NXX+WNP1(3,R)*NYY) / (WNP1(1,R) * a_R)

   !Part 12:
	IF(Dabs(M_L)<1.)Then
     M_Plus = 0.25*(M_L+1.)*(M_L+1.)
	 P_Plus = 0.25*(M_L+1.)*(M_L+1.)*(2.-M_L)
	Else
	 M_Plus = 0.5*(M_L+Dabs(M_L))
	 P_Plus = 0.5*(M_L+Dabs(M_L))/(M_L)
	End If

   !Part 13:
	IF(Dabs(M_R)<1.)Then
     M_Minus =-0.25*(M_R-1.)*(M_R-1.)
	 P_Minus = 0.25*(M_R-1.)*(M_R-1.)*(2.+M_R)
	Else
	 M_Minus = 0.5*(M_R-Dabs(M_R))
	 P_Minus = 0.5*(M_R-Dabs(M_R))/(M_R)
	End If

   !Part 14:
    Mm = M_Plus+M_Minus
	Pm = P(L)*P_Plus + P(R)*P_Minus 

   !Part 15:
	If(Mm>0.)Then
     Ro = WNP1(1,L)       * a_L
	 RU = WNP1(2,L)       * a_L
	 RV = WNP1(3,L)       * a_L
	 RH =(WNP1(4,L)+P(L)) * a_L
	Else
     Ro = WNP1(1,R)       * a_R
	 RU = WNP1(2,R)       * a_R
	 RV = WNP1(3,R)       * a_R
	 RH =(WNP1(4,R)+P(R)) * a_R
	Endif

   !Part 16:
    F1 = ( Mm * Ro          ) * DAA
    F2 = ( Mm * RU + Pm*NXX ) * DAA
    F3 = ( Mm * RV + Pm*NYY ) * DAA
    F4 = ( Mm * RH          ) * DAA

   !Part 17:
    Con(1,L) = Con(1,L) + F1
    Con(2,L) = Con(2,L) + F2
    Con(3,L) = Con(3,L) + F3
    Con(4,L) = Con(4,L) + F4

   !Part 18:
    Con(1,R) = Con(1,R) - F1
    Con(2,R) = Con(2,R) - F2
    Con(3,R) = Con(3,R) - F3
    Con(4,R) = Con(4,R) - F4

 End Do
 
 !pause
!*********************************************************************************************
 End
!###########################################################################################
