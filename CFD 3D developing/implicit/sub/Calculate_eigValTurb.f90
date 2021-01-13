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
 Subroutine Calculate_eigValTurb(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA,eigen)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,WNP1,WB,DA
 Intent (Out  )::eigen

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::R,U,V,W,R_L,RU_L,RV_L,RW_L,R_R,RU_R,RV_R,RW_R,NXX,NYY,NZZ,DAA,eigL,eigR
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,eigen
 Integer,Dimension(1:6,1:Dim)::IDS
!*******************************************************************************************
!Part 1:
 DO I=NF2+1,NF

   !Part 2:
    R = WB(1,I)
    U = WB(2,I)/R
    V = WB(3,I)/R
    W = WB(4,I)/R

   !Part 3:
    eigen(I)=abs( ( U*NX(I) + V*NY(I) + W*NZ(I) )/DA(I) )
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
    eigL=abs( ( RU_L*NXX + RV_L*NYY + RW_L*NZZ ) / (R_L*DAA) )

   !Part 9:
    eigR=abs( ( RU_R*NXX + RV_R*NYY + RW_R*NZZ ) / (R_R*DAA) )

   !Part 10:
    eigen(I)=MAX(eigR,eigL)
 End Do
!*******************************************************************************************
 End
!###########################################################################################

