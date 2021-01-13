!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate the source Terms of Turbulence Model                       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,H Kharinezhad Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F051F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Ke_Source(Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY,St)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,IDS,DW,A,INW,MR,Ce1,Ce2,WNP1,WTNP1,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY
 Intent(Out  )::St

 Integer::Dim,I,II,NC,ME,P1,P2
 Real(8)::K,Epsilon,Rho,MR,Yn,Ce1,Ce2,Pk,Pe,Txx,Txy,Tyy,Lk,Le,fe1,fe2,Tauwall,Ustar,Yplus,Rt
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,DUY,DDUX,DDUY,DDVX,DDVY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
!******************************************************************************************* 
!Part 1:
 Do I=1,NC
           
   !Part 2:
    Rho     = WNP1(1,I)
    k       = WTNP1(1,I)/Rho
    Epsilon = WTNP1(2,I)/Rho
           
   !Part 3:
   !Part 4:     
   !Part 5: 
 
   !Part 6:
    Txx = MR * Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) - Rho*K/1.5
    Txy = MR * Mut(I)*( DDUY(I)+DDVX(I)  )
    Tyy = MR * Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) - Rho*K/1.5

    Pe =  txx*DDUX(I) + txy*(DDUY(I)+DDVX(I)) + tyy*DDVY(I) 
    Pk = min (Pe,10.0*Rho*Epsilon)
 
   !Part 7:  
    St(1,I) = A(I) * ( -Pk + Rho*Epsilon ) 
    St(2,I) = A(I) * ( -Ce1*Epsilon*Pk/k + Ce2*Rho*Epsilon*Epsilon/k )   

 END DO
!*******************************************************************************************
 End
!###########################################################################################

