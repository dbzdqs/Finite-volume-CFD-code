!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:To Initialize all of the Parameters Contribute in Turbulence Model       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar,M.A.Zoljanahi Iran, Tehran, OpenFlows@chmail.ir              //!
!// Doc ID: MC2F060F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_Init3D(Dim,NC,IDS,Minf,Rinf,Wtnp1,Mut,TUinf,Mutinf,&
                               Ginf,Kinf,Oinf,WNP1) 
                               
 Implicit None
!*******************************************************************************************
  Integer::Dim,J,I,II,P1,P2,ME,NC,NFW1,NFW2
 Real(8)::Minf,Rinf,Dmin,Dis,Xj,Yj,Xi,Yi,DX,DY,TUinf,Kinf,Oinf,Mutinf,Ginf
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:3,1:Dim)::Wnt,Wtnp1
 Real(8),Dimension(1:Dim)::Mut
 Real(8),Dimension(1:5,1:Dim)::WNP1
 
 Intent(In )::Dim,NC,Minf,Rinf,IDS
 Intent(Out)::Wtnp1,Mut,Kinf,Oinf,WNP1,Ginf

!*******************************************************************************************
  !Part 1:
  TUinf=0.2
  Mutinf=10.0
  Kinf = (Minf*TUinf/100.0)*(Minf*TUinf/100.0)*1.5 !????????????
  Oinf = 0.09*Kinf/Mutinf*Rinf/Minf  !?????????????????
  Ginf = 1.0
  
  print*, 'Rinf', Rinf ,'TUinf',TUinf, 'Mutinf' ,Mutinf ,'Ginf' ,Ginf
  
  !Part 2:
 Do J=1,NC
    Wtnp1(1,J)=Kinf
    Wtnp1(2,J)=Oinf
    Wtnp1(3,J)=Ginf

    Mut(J) = (Rinf/Minf)*(Wtnp1(1,J)/Wtnp1(2,J)) 
 End Do     

!*******************************************************************************************
 End
!###########################################################################################