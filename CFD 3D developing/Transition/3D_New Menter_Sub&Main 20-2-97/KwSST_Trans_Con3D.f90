!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Convection Terms of 3 eq. Transition Model 3D             //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                           //!
!// Doc ID: MC2F053F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!**********************************************************************************************
 Subroutine KwSST_Trans_Con3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB,Cont)
 Implicit None
!**********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,WNP1,WTNP1,WB,WTB
 Intent(Out  )::Cont

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE,P1,P2
 Real(8)::U,V,W,DX,DY,F1,F2,F3,Q,R,RU,RV,RK,RFi,Rga,RW
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:Dim)::NX,NY,NZ
 Real(8),Dimension(1:3,1:Dim)::WTNP1,WTB,Cont
 Integer,Dimension(1:6,1:Dim)::IDS
!**********************************************************************************************	
!Part 1:
 DO I=1,NC
    Cont(1,I) = 0.0
    Cont(2,I) = 0.0
    Cont(3,I) = 0.0
 End Do
 
!Part 2:
 DO I=NF2+1,NF
 
   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)
    W = WB(4,I)/WB(1,I)

   !Part 5:
    Q  = U*NX(I)+V*NY(I)+W*NZ(I)

   !Part 6:
    F1 = Q*WTB(1,I)
    F2 = Q*WTB(2,I)
    F3 = Q*WTB(3,I)

   !Part 7:
    Cont(1,ME) = Cont(1,ME) + F1
    Cont(2,ME) = Cont(2,ME) + F2
    Cont(3,ME) = Cont(3,ME) + F3

 End Do

!Part 8:
 DO I=NF1+1,NF2

   !Part 9:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 10:  
   	R  = (Wnp1(1,ME)+Wnp1(1,NE))/2.0  
	RU = (Wnp1(2,ME)+Wnp1(2,NE))/2.0
	RV = (Wnp1(3,ME)+Wnp1(3,NE))/2.0
	RW = (Wnp1(4,ME)+Wnp1(4,NE))/2.0
    
    U = RU / R
    V = RV / R
	W = RW / R

   !Part 11:
    Q  = U*NX(I)+V*NY(I)+W*NZ(I)
    
    !Part 12:
    if( Q>=0.0 )then
     RK  = WTNP1(1,ME)
     RFi = WTNP1(2,ME)
     RGa = WTNP1(3,ME)
    else
     RK  = WTNP1(1,NE)
     RFi = WTNP1(2,NE)
     RGa = WTNP1(3,NE)
    end if
    
    F1 = Q * RK
    F2 = Q * RFi
    F3 = Q * RGa

   !Part 13:
    Cont(1,ME) = Cont(1,ME) + F1
    Cont(2,ME) = Cont(2,ME) + F2
    Cont(3,ME) = Cont(3,ME) + F3

    Cont(1,NE) = Cont(1,NE) - F1
    Cont(2,NE) = Cont(2,NE) - F2
    Cont(3,NE) = Cont(3,NE) - F3

 End Do
!*******************************************************************************************
 End
!###########################################################################################
