!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate wall function Model of Turbulence Model                    //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,H.Kharinezhad Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F051F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-ComMErcial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine wallfunc_swfke(Dim,NC,NX,NY,NFW1,NFW2,IDS,DW,A,INW,MR,WNP1,WTNP1,Mu,Mut,&
                                iwf ,x,y,St)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NX,NY,IDS,DW,A,INW,MR,WNP1,Mu,Mut,NFW1,NFW2,x,y
 
 Intent(inOut  )::iwf,WTNP1,St
 Integer::Dim,I,II,NC,ME,P1,P2,tnt,NFW1,NFW2
 Real(8)::K,Epsilon,Rho,MR,Yn,Ce1,Ce2,Pk,Pe,Txx,Txy,Tyy,Lk,Le,fe1,fe2,Tauwall,Ustar,Yplus,Rt ,Cmu,miu
 Real(8)::DX,DY,u,v,Up,Uplus,E,Kapa,c1,c2,dx2,dy2,eps,yplus2,resyp,rm,y11(1:1000),DL
 Real(8)::kplus,eplus,y111(1:100),y10(1:100),y12(1:100),y13(1:100)
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,x,y,NX,NY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,St
 Real(8),Dimension(1:3,1:Dim)::iwf
 
!******************************************************************************************* 
!Part 1: 
 
 tnt =0
 kapa  = 0.41
 E     = 9.8
 eps    =10**(-30)
 c1     = 0.000096
 c2     = 0.2d0
 Ce1=1.35d0
 Ce2=1.8d0
 Cmu=0.09d0

  DO I=NFW1+1,NFW2
      
       !Part 3:
        ME = IDS(1,I)
        P1 = IDS(3,I)
        P2 = IDS(4,I)
            
        !Part 4:
        DX      = ny(I)
        DY      = nx(I)
        Rho     = WNP1(1,ME)
        k       = WTNP1(1,ME)/Rho
        Epsilon = WTNP1(2,ME)/Rho
        miu     = mu(ME)
        Yn      = DW(ME) 
            
        U       = WNP1(2,ME)/WNP1(1,ME)
        V       = WNP1(3,ME)/WNP1(1,ME)
            
        DL      = Dsqrt(DX*DX+DY*DY)  
        Up      = DABS(U*(DX)+V*(DY))/DL
            
        !Part 5:
        Tauwall = Miu*up/yn
        Ustar   = max( Dsqrt(Dabs(tauwall/Rho)),cmu**(0.25)*Dsqrt(k)/(dsqrt(MR)))
        Yplus   = Rho*ustar*Yn/miu/(dsqrt(MR))
                  
        do II =1,20
            
            if (Yplus>11.63) then  
                Uplus   = dlog(E*Yplus)/kapa*(dsqrt(MR))  
                Ustar   = Up/uplus
            else
                Uplus   = Yplus*(dsqrt(MR))
                Ustar   = Up/uplus
            end if

            Yplus   = Rho*ustar*Yn/miu/(dsqrt(MR))
            
        end do

        !Part 8:
        Tauwall = Rho*ustar*UP/Uplus

        iwf(1,I)    = Yplus
        iwf(2,I)    = Tauwall
        iwf(3,I)    = ustar
            
        !Part 9:
        kplus   = ustar*ustar*(MR)
        eplus   = ustar*ustar*ustar*ustar/miu*Rho*dsqrt(MR)
           
        if (yplus>11.63) then
            WTNP1(1,ME) = rho*kplus/DSQRT(cmu)
            WTNP1(2,ME) = miu/ustar/kapa/yn*eplus*(MR)
        else 
            WTNP1(1,ME) = c1*yplus*yplus*rho*kplus
            WTNP1(2,ME) = c2*rho*eplus*dsqrt(MR)
        end if

        !Part 10:
         St(1,ME) = 0d0
         St(2,ME) = 0d0 
            
  end do

!*******************************************************************************************
 End
!###########################################################################################


    