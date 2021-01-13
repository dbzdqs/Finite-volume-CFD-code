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
 Subroutine Ke_Init2(Dim,NC,NF,NFW1,NFW2,IDS,X,Y,Xc,Yc,MR,DW,INW,IWF,WTNP1,Mut,NN,wnp1,mu) 
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,IDS,X,Y,Xc,Yc,MR,NF,wnp1,mu
 Intent(Out  )::DW,INW,WTNP1,Mut,IWF,NN

 Integer::Dim,J,I,II,P1,P2,ME,NC,NFW1,NFW2,NF
 Real(8)::MR,Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY,l0,cbc,u,v,vn,rho
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:2,1:Dim)::WTNP1
   Real(8),Dimension(1:4,1:Dim)::WNP1

 Real(8),Dimension(1:3,1:Dim)::IWF
 Real(8),Dimension(1:Dim)::X,Y,Xc,Yc,DW,Mut,mu
 Integer,Dimension(1:2, 1:Dim)::nn

!*******************************************************************************************	
 
 !Part 1:
 
 l0=0.0005
 cbc=0.003

 !Part 2:

 IWF(:,:)=11.6d0

 DO I=NFW1+1,NFW2
        
    ME      = IDS(1,I)
    IWF(3,I)= Mu(ME)
    
 END DO 
 
 !Part 3: 	

 Do J=1,NC

    Rho        = WNP1(1,J)
    WTNP1(1,J) = Rho*(Mu(J)/Rho/l0)**2d0*MR*MR  
    WTNP1(2,J) = 0.09*WNP1(1,J)*(WTNP1(1,J)/WNP1(1,J))**1.5d0/l0
    
    WTNP1(1,J) = 1.0e-6
    WTNP1(2,J) = 4.5e-7
    Mut(J)     = (WTNP1(1,J)*WTNP1(1,J)/WTNP1(2,J) ) / MR
    
 End Do

 !Part 4:
 Do J=1,NC
    
   !Part 5:
	Xj = Xc(J)
	Yj = Yc(J)
   
   !Part 6:    
    Dmin=1000000.0

   !Part 7:
    Do I=NFW1+1,NFW2
      
	  !Part 8:
       ME = IDS(1,I)
	   P1 = IDS(3,I)
	   P2 = IDS(4,I)
      
	  !Part 9:
       Xi = 0.5*(X(P1) + X(P2))
       Yi = 0.5*(Y(P1) + Y(P2))
      
	  !Part 10:  
	   DX = Xj-Xi
	   DY = Yj-Yi
       Dis = Dsqrt(DX*DX+DY*DY) 
      
	  !Part 11:
	   If(Dis<Dmin)then
	    Dmin = Dis
	    II   = I
	   Endif

	End do
      
   !Part 12: 
    DW(j)  = Dmin
    INW(j) = II

 End do
 
 !Part 13: 	

 call northcellwall(Dim,NF,NC,NFW1,NFW2,IDS,NN,Xc,Yc)
 
!*******************************************************************************************
 End
!###########################################################################################
    
    
!*******************************************************************************************
 Subroutine northcellwall(Dim,NF,NC,NFW1,NFW2,IDS,NN,Xc,Yc)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS,NFW1,NFW2,Xc,Yc
 Intent(Inout)::NN
 
 Integer::Dim,NF,NC,ME,NE,I,J,J1,J2,E,E1,E2,E3,P1_E2,P2_E1,P,ME1,p11,p12
 Integer::NFW1,NFW2,p1,p2
  Real(8)::MR,Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY
 Integer,Dimension(1:4, 1:Dim)::IDS
 Integer,Dimension(1:4, 1:Dim)::CELL_EDGE,Corn
 Integer,Dimension(1:Dim)::NCELL_EDGE
 Integer,Dimension(1:2, 1:Dim)::nn
  Real(8),Dimension(1:Dim)::Xc,Yc
!******************************************************************************************* 
!Part 1:
 Do I=1,NC
    NCELL_EDGE(I) = 0
 End Do
 Corn    = 0

!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1, I)
    NE = IDS(2, I)
	
    !Part 4:
    IF(ME/=0)Then

    !Part 5:
	 NCELL_EDGE(ME) = NCELL_EDGE(ME) + 1
     CELL_EDGE(NCELL_EDGE(ME),ME)=I

    EndIF

	!Part 6:
    IF(NE/=0)Then

     !Part 7:
	 NCELL_EDGE(NE) = NCELL_EDGE(NE) + 1
     CELL_EDGE(NCELL_EDGE(NE),NE)=-I

    EndIF

 End Do


!Part 8:
 Do I=1,NC

    Do J1=1,NCELL_EDGE(I)-1

       E1 = CELL_EDGE(J1,I)

    !Part 9:
    IF( E1>0 )Then
     P2_E1 = IDS(4,E1)
    Else
     P2_E1 = IDS(3,-E1)
	EndIF

	   Do J2=2,NCELL_EDGE(I)
    !Part 10:
    E2 = CELL_EDGE(J2,I)

    !Part 11:
    IF( E2>0 )Then
     P1_E2 = IDS(3,E2)
    Else
     P1_E2 = IDS(4,-E2)
    EndIF


    !Part 12:
    IF( P2_E1==P1_E2 )Then
	 E                 = CELL_EDGE(J1+1,I)
	 CELL_EDGE(J1+1,I) = CELL_EDGE(J2,I)
     CELL_EDGE(J2,I)   = E
    EndIF

       End Do

    End Do

 End Do 



 !Part 13:
 
	  !Part 14:
  Do I=NFW1+1,NFW2
      
    ME = IDS(1,I)
	P1 = IDS(3,I)
	P2 = IDS(4,I)
      
    Do j1=1,NCELL_EDGE(ME)
       j2=abs(CELL_EDGE(j1,ME))
        
       ME1  = IDS(1,j2)
       NE  = IDS(2,j2)
       P11 = IDS(3,j2)
       P12 = IDS(4,j2)

       if((p11/=p1).and.(p11/=p2).and.(p12/=p2).and.(p12/=p1)) then
        NN(1,ME)=NE
        NN(2,ME)=j2
        if (ME==NE) NN(1,ME)=ME1
        end if
        
    end do
    
  end do

  !Part 18:

 Do I=NFW1+1,NFW2
     
	  !Part 21:  
      ME = IDS(1,I)
      Xi = (Xc(nn(1,Me)) - Xc(me))
      Yi = (Yc(nn(1,Me)) - Yc(me))
      
	  !Part 22:  
	   DX = Xi
	   DY = Yi
       Dis = Dsqrt(DX*DX+DY*DY) 

 end do
 
!*******************************************************************************************
 End
!###########################################################################################


