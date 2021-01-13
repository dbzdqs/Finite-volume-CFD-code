
 Subroutine TraSST_Init_V1(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,Dw,INW,R0,Minf,Rinf,Wnt,Wntp1,Mut,TUinf,Mutinf,Init,Wnp1) 
 Implicit None
!*******************************************************************************************
 !Intent(In   )::Dim,NC,R0,Minf,Rinf
 !Intent(Out  )::Wnt,Wntp1,Mut

 Integer::Dim,J,I,II,P1,P2,ME,NC,NFW1,NFW2,Init
 Real(8)::Minf,Rinf,R0,Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY,TUinf,Kinf,Oinf,Mutinf
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::Wnt,Wntp1,Wnp1
 Real(8),Dimension(1:Dim)::X,Y,Xc,Yc,DW,Mut
!*******************************************************************************************	
!Part 1:
 ! TUinf = 0.3
  Kinf = (Minf*TUinf/100.0)*(Minf*TUinf/100.0)*1.5
!  Mutinf = 8.0
  Oinf = Rinf*Kinf/Minf/Mutinf 
  print*, Rinf
  
 Do J=1,NC
    Wntp1(1,J)=Kinf
    Wntp1(2,J)=Oinf
    Wntp1(3,J)=1.0
    if (TUinf>1.3) then
        Wntp1(4,J) = 331.5*((Tuinf - 0.5658)**(-0.671))
    else
        Wntp1(4,J) = (1173.51 - 589.428*Tuinf + 0.2196/Tuinf/Tuinf)
    end if
    Mut(J) = (Rinf/Minf)*(Wntp1(1,J)/Wntp1(2,J)) 
 End Do     

 !Part 2:
 Do J=1,NC
    
   !Part 3:
	Xj = Xc(J)
	Yj = Yc(J)
   
   !Part 4:    
    Dmin=1000000.0

   !Part 5:
    Do I=NFW1+1,NFW2
      
	  !Part 6:
       ME = IDS(1,I)
	   P1 = IDS(3,I)
	   P2 = IDS(4,I)
      
	  !Part 7:
       Xi = 0.5*(X(P1) + X(P2))
       Yi = 0.5*(Y(P1) + Y(P2))
      
	  !Part 8:  
	   DX = Xj-Xi
	   DY = Yj-Yi
       Dis = Dsqrt(DX*DX+DY*DY) 
      
	  !Part 9:
	   If(Dis<Dmin)then
	    Dmin = Dis
	    II   = I
	   Endif

	End do
      
   !Part 10: 
    DW(j)  = Dmin
    INW(j) = II

 End do
 
  IF(Init==1)Then

  Open(4,File='SolutionData.txt')
  
  Do I=1,NC
	Read(4,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I)
 End Do
 
 Do I=1,NC
	Read(4,*) Wntp1(1,I),Wntp1(2,I),Wntp1(3,I),Wntp1(4,I)
 End Do
 
 Do I=1,NC
	Read(4,*) Mut(I)
 End Do      
 
  close(4)

 Endif

 
!*******************************************************************************************
 End