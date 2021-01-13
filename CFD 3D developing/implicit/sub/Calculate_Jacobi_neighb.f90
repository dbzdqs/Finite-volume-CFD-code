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
 Subroutine Calculate_Jacobi_neighb(Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,face,Neib,WNP1,WB,Jaco)
 Implicit None
!*******************************************************************************************
 Intent (In   )::Dim,NF1,NF2,NF,IDS,NX,NY,NZ,DA,GM,face,Neib,WNP1,WB
 Intent (Out  )::Jaco

 Integer::Dim,I,NF,NF1,NF2,ME,NE,k,face,Neib
 Real(8)::NXX,NYY,NZZ,Pm,GM
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ,DA
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:5)::Jaco
 Real(8)::a,h0,ek,Q
 Real(8)::R,U,V,W,H,P
!*******************************************************************************************
!Part 1:
 I=face

!Part 2:
 NE = Neib

!Part 3:
 NXX = NX(I)/DA(I)
 NYY = NY(I)/DA(I)
 NZZ = NZ(I)/DA(I)

!Part 4:
 R = (WNP1(1,NE))
 U = (WNP1(2,NE))/R
 V = (WNP1(3,NE))/R
 W = (WNP1(4,NE))/R
 P = (GM-1)*(WNP1(5,NE)-0.5*R*(U*U+V*V+W*W))
 H = (WNP1(5,NE)+P)/R

!Part 5:
 ek = 0.5*(U**2+V**2+W**2)
 Q  = U*NXX + V*NYY + W*NZZ
 a  = DSQRT((GM-1.0)*(H-ek))
 h0 = (a**2)/(GM-1)+ek

!Part 6:
 Jaco(1,1)=0.0
 Jaco(1,2)=NXX
 Jaco(1,3)=NYY
 Jaco(1,4)=NZZ
 Jaco(1,5)=0.0

 Jaco(2,1)=(GM-1)*ek*NXX-u*Q
 Jaco(2,2)=Q-(GM-2)*u*NXX
 Jaco(2,3)=u*NYY-(GM-1)*v*NXX
 Jaco(2,4)=u*NZZ-(GM-1)*w*NXX
 Jaco(2,5)=(GM-1)*NXX

 Jaco(3,1)=(GM-1)*ek*NYY-v*Q
 Jaco(3,2)=v*NXX-(GM-1)*u*NYY
 Jaco(3,3)=Q-(GM-2)*v*NYY
 Jaco(3,4)=v*NZZ-(GM-1)*w*NYY
 Jaco(3,5)=(GM-1)*NYY

 Jaco(4,1)=(GM-1)*ek*NZZ-w*Q
 Jaco(4,2)=w*NXX-(GM-1)*u*NZZ
 Jaco(4,3)=w*NYY-(GM-1)*V*NZZ
 Jaco(4,4)=Q-(GM-2)*w*NZZ
 Jaco(4,5)=(GM-1)*NZZ

 Jaco(5,1)=((GM-1)*ek-h0)*Q
 Jaco(5,2)=h0*NXX-(GM-1)*u*Q
 Jaco(5,3)=h0*NYY-(GM-1)*v*Q
 Jaco(5,4)=h0*NZZ-(GM-1)*w*Q
 Jaco(5,5)=GM*Q
!*******************************************************************************************
End
!###########################################################################################

