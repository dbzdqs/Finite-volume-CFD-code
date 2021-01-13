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
 Subroutine Calculate_eigMatrixMeanFlow(Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face,eigenMatrix)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NF2,NF,IDS,NX,NY,NZ,DA,GM,WNP1,WB,P,Face
 Intent (Out  )::eigenMatrix

 Integer::Dim,I,NF,NF2,ME,NE,k,kk,Face,cond
 Real(8)::NXX,NYY,NZZ,DAA,Pm,GM
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA,P
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:5)::eigenMatrix,RR,LL,Landa,CC
 Real(8),Dimension(1:5)::Dumy
 Real(8)::a,h0,ek,Q,R_RT
 Real(8)::R_L,U_L,V_L,W_L,H_L,A_L,R_R,U_R,V_R,W_R,H_R,A_R,R,U,V,W,H,P_L,P_R,Q_L,Q_R
!*******************************************************************************************

!Part 1:
 I=Face

!Part 2:
 ME = IDS(1,I)
 NE = IDS(2,I)

!Part 3:
 DAA = DA(I)
 NXX = NX(I)/DAA
 NYY = NY(I)/DAA
 NZZ = NZ(I)/DAA

!Part 4:
 If(NF2<Face .and. Face<NF+1) then

 !Part 5:
  R  = WB(1,I)
  U  = WB(2,I)/WB(1,I)
  V  = WB(3,I)/WB(1,I)
  W  = WB(4,I)/WB(1,I)
  Pm = WB(6,I)
  H  = (WB(5,I)+Pm)/R

!Part 6:
 Else

 !Part 7:
  R_L = WNP1(1,ME)
  U_L = (WNP1(2,ME))/R_L
  V_L = (WNP1(3,ME))/R_L
  W_L = (WNP1(4,ME))/R_L

 !Part 8:
  R_R = WNP1(1,NE)
  U_R = (WNP1(2,NE))/R_R
  V_R = (WNP1(3,NE))/R_R
  W_R = (WNP1(4,NE))/R_R

 !Part 9:
  P_L = (GM-1)*(WNP1(5,ME)-0.5*R_L*(U_L*U_L+V_L*V_L+W_L*W_L))
  P_R = (GM-1)*(WNP1(5,NE)-0.5*R_R*(U_R*U_R+V_R*V_R+W_R*W_R))

 !Part 10:
  H_L = (WNP1(5,ME)+P_L)/R_L
  H_R = (WNP1(5,NE)+P_R)/R_R

 !Part 11:
  A_L  = DSQRT(GM*P_L/R_L)
  A_R  = DSQRT(GM*P_R/R_R)

 !Part 12:
  R_RT = DSQRT(R_R/R_L)

 !Part 13:
  R = DSQRT(R_L*R_R)
  U = (U_L+U_R*R_RT)/(1.0+R_RT)
  V = (V_L+V_R*R_RT)/(1.0+R_RT)
  W = (W_L+W_R*R_RT)/(1.0+R_RT)
  H = (H_L+H_R*R_RT)/(1.0+R_RT)
 End if

!Part 14:
 ek = 0.5*(U**2+V**2+W**2)
 Q  = U*NXX + V*NYY + W*NZZ
 a  = DSQRT((GM-1.0)*(H-ek))
 h0 = (a**2)/(GM-1)+ek

!Part 15:
 Landa(:,:)=0.0

!Part 16:
 If(NF2<Face .and. Face<NF+1) then
  Landa(1,1)=abs(Q)+a
  Landa(2,2)=abs(Q)+a
  Landa(3,3)=abs(Q)+a
  Landa(4,4)=abs(Q)+a
  Landa(5,5)=abs(Q)+a

!Part 17:
 Else
  Landa(1,1)=Q-a
  Landa(2,2)=Q
  Landa(3,3)=Q+a
  Landa(4,4)=Q
  Landa(5,5)=Q
 End if

!Part 18:
!Landa(:,:)=0.5*(Landa(:,:)+abs(Landa(:,:)))
Landa(:,:)=abs(Landa(:,:))

!Part 19:
 If(abs(NXX)>0.001)then
  cond=1

 !Part 20:
  RR(1,1)=1.0
  RR(1,2)=1.0
  RR(1,3)=1.0
  RR(1,4)=0.0
  RR(1,5)=0.0

  RR(2,1)=u-a*NXX
  RR(2,2)=u
  RR(2,3)=u+a*NXX
  RR(2,4)=NYY
  RR(2,5)=-NZZ

  RR(3,1)=v-a*NYY
  RR(3,2)=v
  RR(3,3)=v+a*NYY
  RR(3,4)=-NXX
  RR(3,5)=0.0

  RR(4,1)=w-a*NZZ
  RR(4,2)=w
  RR(4,3)=w+a*NZZ
  RR(4,4)=0
  RR(4,5)=NXX

  RR(5,1)=h0-a*Q
  RR(5,2)=ek
  RR(5,3)=h0+a*Q
  RR(5,4)=u*NYY-v*NXX
  RR(5,5)=w*NXX-u*NZZ

 !Part 21:
  LL(1,1)=((GM-1)*ek+a*Q)/(2*a**2)
  LL(1,2)=(-(GM-1)*u-a*NXX)/(2*a**2)
  LL(1,3)=(-(GM-1)*v-a*NYY)/(2*a**2)
  LL(1,4)=(-(GM-1)*w-a*NZZ)/(2*a**2)
  LL(1,5)=((GM-1))/(2*a**2)

  LL(2,1)=(a**2-(GM-1)*ek)/(a**2)
  LL(2,2)=((GM-1)*u)/(a**2)
  LL(2,3)=((GM-1)*v)/(a**2)
  LL(2,4)=((GM-1)*w)/(a**2)
  LL(2,5)=(-(GM-1))/(a**2)

  LL(3,1)=((GM-1)*ek-a*Q)/(2*a**2)
  LL(3,2)=(-(GM-1)*u+a*NXX)/(2*a**2)
  LL(3,3)=(-(GM-1)*v+a*NYY)/(2*a**2)
  LL(3,4)=(-(GM-1)*w+a*NZZ)/(2*a**2)
  LL(3,5)=((GM-1))/(2*a**2)

  LL(4,1)=(v-Q*NYY)/NXX
  LL(4,2)=NYY
  LL(4,3)=(NYY**2-1)/NXX
  LL(4,4)=(NYY*NZZ)/NXX
  LL(4,5)=0.0

  LL(5,1)=(Q*NZZ-w)/NXX
  LL(5,2)=-NZZ
  LL(5,3)=-(NYY*NZZ)/NXX
  LL(5,4)=(1-NZZ**2)/NXX
  LL(5,5)=0.0

!Part 22:
 Else If(abs(NZZ)>0.001)then
  cond=3

 !Part 23:
  RR(1,1)=1.0
  RR(1,2)=1.0
  RR(1,3)=1.0
  RR(1,4)=0.0
  RR(1,5)=0.0

  RR(2,1)=u-a*NXX
  RR(2,2)=u
  RR(2,3)=u+a*NXX
  RR(2,4)=-NZZ
  RR(2,5)=0.0

  RR(3,1)=v-a*NYY
  RR(3,2)=v
  RR(3,3)=v+a*NYY
  RR(3,4)=0.0
  RR(3,5)=NZZ

  RR(4,1)=w-a*NZZ
  RR(4,2)=w
  RR(4,3)=w+a*NZZ
  RR(4,4)=NXX
  RR(4,5)=-NYY

  RR(5,1)=h0-a*Q
  RR(5,2)=ek
  RR(5,3)=h0+a*Q
  RR(5,4)=w*NXX-u*NZZ
  RR(5,5)=v*NZZ-w*NYY

 !Part 24:
  LL(1,1)=((GM-1)*ek+a*Q)/(2*a**2)
  LL(1,2)=(-(GM-1)*u-a*NXX)/(2*a**2)
  LL(1,3)=(-(GM-1)*v-a*NYY)/(2*a**2)
  LL(1,4)=(-(GM-1)*w-a*NZZ)/(2*a**2)
  LL(1,5)=((GM-1))/(2*a**2)

  LL(2,1)=(a**2-(GM-1)*ek)/(a**2)
  LL(2,2)=((GM-1)*u)/(a**2)
  LL(2,3)=((GM-1)*v)/(a**2)
  LL(2,4)=((GM-1)*w)/(a**2)
  LL(2,5)=(-(GM-1))/(a**2)

  LL(3,1)=((GM-1)*ek-a*Q)/(2*a**2)
  LL(3,2)=(-(GM-1)*u+a*NXX)/(2*a**2)
  LL(3,3)=(-(GM-1)*v+a*NYY)/(2*a**2)
  LL(3,4)=(-(GM-1)*w+a*NZZ)/(2*a**2)
  LL(3,5)=((GM-1))/(2*a**2)

  LL(4,1)=(u-Q*NXX)/NZZ
  LL(4,2)=(NXX**2-1)/NZZ
  LL(4,3)=(NXX*NYY)/NZZ
  LL(4,4)=NXX
  LL(4,5)=0.0

  LL(5,1)=(Q*NYY-v)/NZZ
  LL(5,2)=(-NXX*NYY)/NZZ
  LL(5,3)=(1-NYY**2)/NZZ
  LL(5,4)=-NYY
  LL(5,5)=0.0

!Part 25:
 Else
  cond=2

 !Part 26:
  RR(1,1)=1.0
  RR(1,2)=1.0
  RR(1,3)=1.0
  RR(1,4)=0.0
  RR(1,5)=0.0

  RR(2,1)=u-a*NXX
  RR(2,2)=u
  RR(2,3)=u+a*NXX
  RR(2,4)=NYY
  RR(2,5)=0.0

  RR(3,1)=v-a*NYY
  RR(3,2)=v
  RR(3,3)=v+a*NYY
  RR(3,4)=-NXX
  RR(3,5)=NZZ

  RR(4,1)=w-a*NZZ
  RR(4,2)=w
  RR(4,3)=w+a*NZZ
  RR(4,4)=0
  RR(4,5)=-NYY

  RR(5,1)=h0-a*Q
  RR(5,2)=ek
  RR(5,3)=h0+a*Q
  RR(5,4)=u*NYY-v*NXX
  RR(5,5)=v*NZZ-w*NYY

 !Part 27:
  LL(1,1)=((GM-1)*ek+a*Q)/(2*a**2)
  LL(1,2)=(-(GM-1)*u-a*NXX)/(2*a**2)
  LL(1,3)=(-(GM-1)*v-a*NYY)/(2*a**2)
  LL(1,4)=(-(GM-1)*w-a*NZZ)/(2*a**2)
  LL(1,5)=((GM-1))/(2*a**2)

  LL(2,1)=(a**2-(GM-1)*ek)/(a**2)
  LL(2,2)=((GM-1)*u)/(a**2)
  LL(2,3)=((GM-1)*v)/(a**2)
  LL(2,4)=((GM-1)*w)/(a**2)
  LL(2,5)=(-(GM-1))/(a**2)

  LL(3,1)=((GM-1)*ek-a*Q)/(2*a**2)
  LL(3,2)=(-(GM-1)*u+a*NXX)/(2*a**2)
  LL(3,3)=(-(GM-1)*v+a*NYY)/(2*a**2)
  LL(3,4)=(-(GM-1)*w+a*NZZ)/(2*a**2)
  LL(3,5)=((GM-1))/(2*a**2)

  LL(4,1)=(Q*NXX-u)/NYY
  LL(4,2)=(1-NXX**2)/NYY
  LL(4,3)=-NXX
  LL(4,4)=(-NXX*NZZ)/NYY
  LL(4,5)=0.0

  LL(5,1)=(w-Q*NZZ)/NYY
  LL(5,2)=(NXX*NZZ)/NYY
  LL(5,3)=NZZ
  LL(5,4)=(NZZ**2-1)/NYY
  LL(5,5)=0.0
 End If

!Part 28:
 Call Multiply(RR,Landa,CC)

!Part 29:
 Call Multiply(CC,LL,eigenMatrix)

!*******************************************************************************************
End
!###########################################################################################

