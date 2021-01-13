Subroutine TraSST_BC_V1(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,Wbt,Wnp1,Wntp1,IDS,Minf,Rinf,MR,Mu,TUinf,Mutinf)
 Implicit None
!**********************************************************************************************
 !Intent (In   )::Dim,NFW,NFF
 !Intent (Out  )::Wbt

 Integer::Dim,I,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,ME,P1,P2
 Real(8)::MR,U,V,Q,Minf,Rinf,TUinf,Kinf,Oinf,Mutinf
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Mu,DW,NX,NY
 Real(8),Dimension(1:4,1:Dim)::Wbt,Wntp1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:4,1:Dim)::Wnp1
!**********************************************************************************************	
  !Part 1:
 ! TUinf = 0.3
  Kinf = (Minf*TUinf/100.0)*(Minf*TUinf/100.0)*1.5
  !Mutinf = 8.0
  Oinf = Rinf*Kinf/Minf/Mutinf
  
  DO I=NFI1+1,NFI2
 	Wbt(1,I) = Kinf
    Wbt(2,I) = Oinf 
    Wbt(3,I) = 1.0
    if (TUinf>1.3) then
        Wbt(4,I) = 331.5*((Tuinf - 0.5658)**(-0.671))
    else
        Wbt(4,I) = (1173.51 - 589.428*Tuinf + 0.2196/Tuinf/Tuinf)
    end if
    
  END do
  
  !Part 2:  
  DO I=NFO1+1,NFO2
    ME  = IDS(1,I)
 	Wbt(1,I) = Wntp1(1,ME)
    Wbt(2,I) = Wntp1(2,ME)
    Wbt(3,I) = Wntp1(3,ME)
    Wbt(4,I) = Wntp1(4,ME)
  END do
  
  !Part 3:
  DO I=NFW1+1,NFW2
    ME  = IDS(1,I)
 	Wbt(1,I) = 0.0
    Wbt(2,I) = Wb(1,I)*(MR)*(60.0*Mu(ME)/ (Wnp1(1,ME)*0.075*DW(ME)*DW(ME)) ) 
    Wbt(3,I) = Wntp1(3,ME)
    Wbt(4,I) = Wntp1(4,ME)
  END do
  
  !Part 4:
  DO I=NFS1+1,NFS2
    ME  = IDS(1,I)
 	Wbt(1,I) = Wntp1(1,ME)
    Wbt(2,I) = Wntp1(2,ME)
    Wbt(3,I) = Wntp1(3,ME)
    Wbt(4,I) = Wntp1(4,ME)
  END do
  
  !part 5:
  DO I=NFF1+1,NFF2
    ME  = IDS(1,I)
    
    U = Wb(2,I)/Wb(1,I)
    V = Wb(3,I)/Wb(1,I)
    
    Q  = U*NX(I)+V*NY(I)
        
    if (Q<=0.0) then
 	Wbt(1,I) = Kinf
    Wbt(2,I) = Oinf
    Wbt(3,I) = 1.0
    if (TUinf>1.3) then
        Wbt(4,I) = 331.5*((Tuinf - 0.5658)**-0.671)
    else
        Wbt(4,I) = (1173.51 - 589.428*Tuinf + 0.2196/Tuinf/Tuinf)
    end if
    else
 	Wbt(1,I) = Wntp1(1,ME)
    Wbt(2,I) = Wntp1(2,ME) 
    Wbt(3,I) = Wntp1(3,ME)
    Wbt(4,I) = Wntp1(4,ME)    
    end if
  END do

  

  

 !****************************************************************************************

!**********************************************************************************************
 End