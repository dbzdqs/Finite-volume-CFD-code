!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: !
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F000F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Calculate_eigValMeanFlow(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,X,Y,Z,xc,yc,zc,MR,eigen)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,X,Y,Z,MR
 Intent (Out  )::eigen

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::U,V,W,NXX,NYY,NZZ,F1,F2,F3,F4,Q,Pm,R,RU,RV,RW,RE,GM,MumL,DL,MumLt,LamP,LamN,PrL,PrT,MR,eigR,eigL
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P,eigen,Mu,Mut,xc,yc,zc
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,Z
!*******************************************************************************************
!Part 1:
 DO I=NF2+1,NF

   !Part 2:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 3:
    R = WB(1,I)
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)
    Pm = WB(6,I)

   !Part 4:
    Q  = U*NX(I)/DA(I) + V*NY(I)/DA(I) + W*NZ(I)/DA(I)

   !Part 5:
    eigen(I)=abs(Q) + sqrt(GM*Pm/WB(1,I))+( Mu(ME)/PrL+Mut(ME)/PrT )*DA(I)*DA(I)/WB(1,I)
 End Do

!Part 6:
 DO I=NF1+1,NF2

   !Part 7:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 8:
    NXX = NX(I)/DA(I)
    NYY = NY(I)/DA(I)
    NZZ = NZ(I)/DA(I)

   !Part 9:
    DL=SQRT((NXX*(xc(ME)-xc(NE)))**2+(NYY*(yc(ME)-yc(NE)))**2+(NZZ*(Zc(ME)-Zc(NE)))**2)

   !Part 10:
    R  = WNP1(1,ME)
	RU = WNP1(2,ME)
	RV = WNP1(3,ME)
	RW = WNP1(4,ME)
	RE = WNP1(5,ME)
    Pm = P(ME)

   !Part 11:
    U = RU / R
    V = RV / R
    W = RW / R
    MumL = Mu(ME)
    MumLt = Mut(ME)

   !Part 12:
    Q = U*NXX + V*NYY + W*NZZ

   !Part 13:
    LamP=abs((Q + sqrt(GM*Pm/R))) + abs(2*MR*(MumL+MumLt)/((DL*R)))  !+ DA(I)*AMAX1(4.0/3.0,GM)*MR*( MumL/PrL+MumLt/PrT ) / (DL*R)
    LamN=abs((Q - sqrt(GM*Pm/R))) + abs(2*MR*(MumL+MumLt)/((DL*R)))  !+ DA(I)*AMAX1(4.0/3.0,GM)*MR*( MumL/PrL+MumLt/PrT ) / (DL*R)

    eigL=MAX(LamP,LamN)

   !Part 14:
    R  = WNP1(1,NE)
	RU = WNP1(2,NE)
	RV = WNP1(3,NE)
	RW = WNP1(4,NE)
	RE = WNP1(5,NE)
    Pm = P(NE)

   !Part 15:
    U = RU / R
    V = RV / R
    W = RW / R
    MumL = Mu(NE)
    MumLt = Mut(NE)

   !Part 16:
    Q = U*NXX + V*NYY + W*NZZ

   !Part 17:
    LamP=abs((Q + sqrt(GM*Pm/R))) + abs(2*MR*(MumL+MumLt)/((DL*R)))  !+ DA(I)*AMAX1(4.0/3.0,GM)*MR*( MumL/PrL+MumLt/PrT ) / (DL*R)
    LamN=abs((Q - sqrt(GM*Pm/R))) + abs(2*MR*(MumL+MumLt)/((DL*R)))  !+ DA(I)*AMAX1(4.0/3.0,GM)*MR*( MumL/PrL+MumLt/PrT ) / (DL*R)
    eigR=MAX(LamP,LamN)

   !Part 18:
    eigen(I)=MAX(eigR,eigL)

 End Do
!*******************************************************************************************
 End
!###########################################################################################

