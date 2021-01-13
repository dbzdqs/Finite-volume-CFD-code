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
 Subroutine Increment(GM,NX,NY,NZ,W1,W2,W3,W4,W5,DW1,DW2,DW3,DW4,DW5,Con_j)
 Implicit None
!*******************************************************************************************
 Intent(In   )::GM,NX,NY,NZ,W1,W2,W3,W4,W5,DW1,DW2,DW3,DW4,DW5
 Intent(Out  )::Con_j

 Integer::i
 Real(8)::GM,NX,NY,NZ,Pm,U,V,W,Q,W1,W2,W3,W4,W5,DW1,DW2,DW3,DW4,DW5,DDW1,DDW2,DDW3,DDW4,DDW5,P
 Real(8),Dimension(1:5)::Con_j
!*******************************************************************************************
!Part 1:
 U = W2/W1
 V = W3/W1
 W = W4/W1

 Pm = (GM-1)*(W5-0.5*W1*(U*U+V*V+W*W))
 
!Part 2:
 Q  = U*NX + V*NY + W*NZ

!Part 3:
 Con_j(1) = Q*W1
 Con_j(2) = Q*W2 + Pm*NX
 Con_j(3) = Q*W3 + Pm*NY
 Con_j(4) = Q*W4 + Pm*NZ
 con_j(5) = Q*W5 + Pm*Q

!Part 4:
 DDW1=W1+DW1
 DDW2=W2+DW2
 DDW3=W3+DW3
 DDW4=W4+DW4
 DDW5=W5+DW5

!Part 5:
 U = DDW2/DDW1
 V = DDW3/DDW1
 W = DDW4/DDW1
 
 P = (GM-1)*(DDW5-0.5*DDW1*(U*U+V*V+W*W))

!Part 6:
 Q  = U*NX + V*NY + W*NZ

!Part 7:
 Con_j(1) = Q*DDW1            - Con_j(1)
 Con_j(2) = Q*DDW2 + P*NX     - Con_j(2)
 Con_j(3) = Q*DDW3 + P*NY     - Con_j(3)
 Con_j(4) = Q*DDW4 + P*NZ     - Con_j(4)
 Con_j(5) = Q*DDW5 + P*Q      - Con_j(5)
!*******************************************************************************************
 End
!###########################################################################################


