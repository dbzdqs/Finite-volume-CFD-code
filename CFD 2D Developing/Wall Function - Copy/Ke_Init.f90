!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Initialize all of the Parameters Contribute in Turbulence Model       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,H Kharinezhad Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F049F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Ke_Init(Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,MR,DW,INW,WTNP1,Mut) 
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,MR
 Intent(Out  )::DW,INW,WTNP1,Mut

 Integer::Dim,J,I,II,P1,P2,ME,NC,NFW1,NFW2
 Real(8)::MR,Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::X,Y,Xc,Yc,DW,Mut
!*******************************************************************************************	
!Part 1: 	
 Do J=1,NC
    WTNP1(1,J) = 1.0e-6
    WTNP1(2,J) = 4.5e-7
    Mut(J)     = (WTNP1(1,J)*WTNP1(1,J)/WTNP1(2,J) ) / MR
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
!*******************************************************************************************
 End
!###########################################################################################