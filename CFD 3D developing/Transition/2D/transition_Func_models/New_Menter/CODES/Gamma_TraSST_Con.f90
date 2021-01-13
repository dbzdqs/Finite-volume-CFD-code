!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// calculation of Convection term for transition model                                  //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
!**********************************************************************************************
 Subroutine Gamma_TraSST_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,Wntp1,Wb,Wbt,Cont)
 Implicit None
!**********************************************************************************************

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE,P1,P2
 Real(8)::U,V,DX,DY,F1,F2,F3,F4,Q,R,RU,RV,RK,ROmeg,RG,RR
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:Dim)::NX,NY
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Wbt,Cont
 Integer,Dimension(1:4,1:Dim)::IDS
!**********************************************************************************************	
!Part 1:
 DO I=1,NC
    Cont(3,I) = 0.0
 End Do

!Part 2:
 DO I=NFW1+1,NFW2
     
    !Part 3:
    ME = IDS(1,I)

   !Part 4:
    Cont(3,ME) = 0.0
 
 End Do

!Part 5:
 DO I=NFW2+1,NF
 
   !Part 6:
    ME = IDS(1,I)

   !Part 7:
    U = Wb(2,I)/Wb(1,I)
    V = Wb(3,I)/Wb(1,I)

   !Part 8:
    Q  = U*NX(I)+V*NY(I)

   !Part 9:
    F3 = Q*Wbt(3,I)

   !Part 10:
    Cont(3,ME) = Cont(3,ME) + F3
 End Do

 
!Part 11:
 DO I=NF1+1,NF2

   !Part 12:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 13:
    R  = (Wnp1(1,ME)+Wnp1(1,NE))/2.0  
	RU = (Wnp1(2,ME)+Wnp1(2,NE))/2.0
	RV = (Wnp1(3,ME)+Wnp1(3,NE))/2.0
       
    U = RU / R
    V = RV / R

   !Part 14:
    Q  = U*NX(I)+V*NY(I)
    
    !Part 15:
    if (Q>=0.0) then
    
    RG=Wntp1(3,ME)
    else
    RG=Wntp1(3,NE)
    end if
    
    F1 = Q * RK
    F2 = Q * ROmeg
    F3 = Q * RG


   !Part 16:
    Cont(3,ME) = Cont(3,ME) + F3

    Cont(3,NE) = Cont(3,NE) - F3
 End Do



!**********************************************************************************************
 End
!##############################################################################################
