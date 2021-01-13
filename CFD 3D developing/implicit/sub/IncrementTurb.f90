!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:  !
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F009F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine IncrementTurb(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA,JacobiL,JacobiR)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA
 Intent(Out  )::JacobiL,JacobiR

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::U,V,W,R_L,RU_L,RV_L,RW_L,R_R,RU_R,RV_R,RW_R,NXX,NYY,NZZ,DAA
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,JacobiL,JacobiR
 Integer,Dimension(1:6,1:Dim)::IDS
!*******************************************************************************************
!Part 1:
 DO I=NF2+1,NF

   !Part 2:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)

   !Part 3:
    JacobiL(I)=( U*NX(I) + V*NY(I) + W*NZ(I) )/DA(I)
    JacobiR(I)=0.0
 End Do

!Part 4:
 DO I=NF1+1,NF2

   !Part 5:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 6:
    DAA = DA(I)
    NXX = NX(I)
    NYY = NY(I)
    NZZ = NZ(I)

   !Part 7:
    R_L  = WNP1(1,ME)
	RU_L = WNP1(2,ME)
	RV_L = WNP1(3,ME)
	RW_L = WNP1(4,ME)

    R_R  = WNP1(1,NE)
	RU_R = WNP1(2,NE)
	RV_R = WNP1(3,NE)
	RW_R = WNP1(4,NE)

   !Part 8:
    JacobiL(I)=( RU_L*NXX + RV_L*NYY + RW_L*NZZ ) / (R_L*DAA)

   !Part 9:
    JacobiR(I)=( RU_R*NXX + RV_R*NYY + RW_R*NZZ ) / (R_R*DAA)
 End Do
!*******************************************************************************************
 End
!###########################################################################################

