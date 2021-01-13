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
 Subroutine Jaco_viscous(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                            &,MR,Cell,Face,Jacobv)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Mu,Mut,PrL,PrT,xc,yc,zc,Vol&
                            &,MR,Cell,Face
 Intent (Out  )::Jacobv

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE,cell,Cell_Neib
 Real(8)::U,V,W,NXX,NYY,NZZ,F1,F2,F3,F4,Q,Pm,R,RU,RV,RW,RE,GM,MumL,DL,MumLt,LamP,LamN,PrL,PrT,MR,eigR,eigL,Volume,SOR
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P,eigen,Mu,Mut,xc,yc,zc,Vol
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:5)::Jacobv
 Integer::j,k,kk,Neib,Face
!*******************************************************************************************
!Part 1:
 Jacobv(:,:)=0.0

!Part 2:
 ME=IDS(1,Face)
 NE=IDS(2,Face)

!Part 3:
 If(NF2<Face .and. Face<NF+1) then

 !Part 4:
  I  = Face
  R  = WB(1,I)
  U  = WB(2,I)/WB(1,I)
  V  = WB(3,I)/WB(1,I)
  W  = WB(4,I)/WB(1,I)
  Pm = WB(6,I)

 !Part 5:
  MumL  = Mu (ME)
  MumLt = Mut(ME)

!Part 6:
 Else

 !Part 7:
  R  = WNP1(1,NE)
  RU = WNP1(2,NE)
  RV = WNP1(3,NE)
  RW = WNP1(4,NE)
  RE = WNP1(5,NE)

 !Part 8:
  U = RU/R
  V = RV/R
  W = RW/R
  Pm = (GM-1)*(RE-0.5*R*(U*U+V*V+W*W))

 !Part 9:
  R  = 0.5*(WNP1(1,ME)+R)
  RU = 0.5*(WNP1(2,ME)+RU)
  RV = 0.5*(WNP1(3,ME)+RV)
  RW = 0.5*(WNP1(4,ME)+RW)
  RE = 0.5*(WNP1(5,ME)+RE)
  Pm = 0.5*(P(ME)+Pm)

 !Part 10:
  U = RU / R
  V = RV / R
  W = RW / R

 !Part 11:
  MumL  = 0.5*(Mu(ME) +Mu (NE))
  MumLt = 0.5*(Mut(ME)+Mut(NE))
 End if

!Part 12:
 Volume=Vol(Cell)

!Part 13:
 NXX = NX(Face)/DA(Face)
 NYY = NY(Face)/DA(Face)
 NZZ = NZ(Face)/DA(Face)

!Part 14:
 if(NE/=0) DL=SQRT((NXX*(xc(ME)-xc(NE)))**2+(NYY*(yc(ME)-yc(NE)))**2+(NZZ*(Zc(ME)-Zc(NE)))**2)
 DL=(MR*(MumL+MumLt)*DA(Face))/DL
 if(NE==0) DL=0.0
 
!Part 15:
 Jacobv(2,2)=Jacobv(2,2)+( (1.0/3.0)*NXX*NXX  + 1)*DL
 Jacobv(2,3)=Jacobv(2,3)+( (1.0/3.0)*NXX*NYY    )*DL
 Jacobv(2,4)=Jacobv(2,4)+( (1.0/3.0)*NXX*NZZ    )*DL

 Jacobv(3,2)=Jacobv(3,2)+( (1.0/3.0)*NXX*NYY    )*DL
 Jacobv(3,3)=Jacobv(3,3)+( (1.0/3.0)*NYY**2   +1)*DL
 Jacobv(3,3)=Jacobv(3,4)+( (1.0/3.0)*NYY*NZZ    )*DL

 Jacobv(4,2)=Jacobv(4,2)+( (1.0/3.0)*NXX*NZZ    )*DL
 Jacobv(4,3)=Jacobv(4,3)+( (1.0/3.0)*NYY*NZZ    )*DL
 Jacobv(4,4)=Jacobv(4,4)+( (1.0/3.0)*NZZ**2   +1)*DL

 Jacobv(5,1)=Jacobv(5,1)+((Pm*GM)/(R*R*PrT*(GM-1)*(MumL+MumLt)))*DL
 Jacobv(5,2)=Jacobv(5,2)+((1.0/3.0)*NXX*Volume+u)*DL
 Jacobv(5,3)=Jacobv(5,3)+((1.0/3.0)*NYY*Volume+v)*DL
 Jacobv(5,4)=Jacobv(5,4)+((1.0/3.0)*NZZ*Volume+w)*DL
 Jacobv(5,5)=Jacobv(5,5)+((GM)/(R*PrT*(GM-1)*(MumL+MumLt)))*DL
!*******************************************************************************************
 End
!###########################################################################################

