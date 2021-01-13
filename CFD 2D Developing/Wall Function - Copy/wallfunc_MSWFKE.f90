!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: To Calculate the source Terms of Turbulence Model                       //!
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
 Subroutine wallfunc_MSWFKe(Dim,NC,NX,NY,NFW1,NFW2,IDS,DW,A,INW,MR,WNP1,WTNP1,Mu,Mut,&
                                IWF ,x,y,NN,st)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NX,NY,IDS,DW,A,INW,MR,WNP1,Mu,Mut,NFW1,NFW2,x,y,NN
 
 Intent(inOut  )::IWF,WTNP1,st
 
 
 Integer::Dim,I,II,NC,ME,P1,P2,tnt,NFW1,NFW2,PN1,PN2,MEN,FNN
 Real(8)::K,Epsilon,Rho,MR,yp,Ce1,Ce2,Pk,Pe,Txx,Txy,Tyy,Lk,Le,fe1,fe2,Tauwall,Ustar,Yplus,Rt ,Cmu,miu
 Real(8)::DX,DY,u,v,Up,Uplus,E,Kapa,c1,c2,dx2,dy2,eps,yplus2,resyp,rm,y11(1:1000),DL,yv,yn
 Real(8)::kplus,eplus,y111(1:100),y10(1:100),y12(1:100),y13(1:100),xi,xj,yj,yi,kN
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut,x,y,NX,NY
 Real(8),Dimension(1:2,1:Dim)::WTNP1,st
 Real(8),Dimension(1:3,1:Dim)::IWF
  Integer,Dimension(1:2, 1:Dim)::nn

!******************************************************************************************* 
!Part 1: 
 
 tnt    = 0
 kapa   = 0.41
 E      = 9.8
 eps    = 10**(-30)
 c1     = 0.000096d0
 c2     = 0.2d0
 Ce1    = 1.35d0
 Ce2    = 1.8d0
 Cmu    = 0.09d0
 !Part 2: 
 
 DO I=NFW1+1,NFW2
  
    ME      = IDS(1,I)
    P1      = IDS(3,I)
    P2      = IDS(4,I)
    FNN     = NN(2,ME)
    MEN     = NN(1,ME)
    PN1     = IDS(3,FNN)
    PN2     = IDS(4,FNN)

    !Part 4:
    DX      = ny(I)
    DY      = nx(I)
    Rho     = WNP1(1,ME)
    k       = WTNP1(1,ME)/Rho
    Epsilon = WTNP1(2,ME)/Rho
    miu     = mu(ME)
    yp      = DW(ME) 
    kN      = WTNP1(1,MEN)/WNP1(1,MEN)
    U       = WNP1(2,ME)/WNP1(1,ME)
    V       = WNP1(3,ME)/WNP1(1,ME)
    DL      = Dsqrt(DX*DX+DY*DY)  
    Up      = DABS(U*(DX)+V*(DY))/DL
            
    !Part 5:
    Xi      = 0.5*(X(P1) + X(P2))
    Yi      = 0.5*(Y(P1) + Y(P2))
    Xj      = 0.5*(X(Pn1) + X(Pn2))
    Yj      = 0.5*(Y(Pn1) + Y(Pn2))
	DX      = Xj-Xi
	DY      = Yj-Yi
    yn      = Dsqrt(DX*DX+DY*DY)       
       
    !Part 6:
    Tauwall = Miu*up/yp
    Ustar   = max( Dsqrt(Dabs(tauwall/Rho)),cmu**(0.25)*Dsqrt(k))
    Yplus   = Rho*ustar*yp/miu/(dsqrt(MR))
            
    !Part 7:
    do II =1,20
            
        if (Yplus>11.63) then    
            Uplus   = dlog(E*Yplus)/kapa*(dsqrt(MR))  
            Ustar   = Up/uplus   
        else         
            Uplus   = Yplus*(dsqrt(MR))
            Ustar   = Up/uplus  
        end if
        
        Yplus   = Rho*ustar*yp/miu/(dsqrt(MR))
        Tauwall = Rho*ustar*ustar  
        
    end do

    !Part 8:
    IWF(1,I) = Yplus
    IWF(2,I) = Tauwall

    !Part 9:
    kplus   = ustar*ustar*(MR)
    eplus   = ustar*ustar*ustar*ustar/miu*Rho*dsqrt(MR)
    yv      = 11.63d0*miu/Rho/dsqrt(k)*((MR))
            
    if (yplus>11.63) then
        WTNP1(1,ME) = rho*kplus/DSQRT(cmu)     
        k           =   WTNP1(1,ME)/Rho
        WTNP1(2,ME) = miu/ustar/kapa/yn*eplus*(MR)
        IWF(3,I)    = IWF(2,I)*yp/up
    else          
        WTNP1(1,ME) = c1*yplus*yplus*rho*kplus
        k           =   WTNP1(1,ME)/Rho
        WTNP1(2,ME) = 2d0*miu*kN/yn/yv*MR+rho*kN**1.5d0/2.44*log(yn/yv)/yn
        WTNP1(2,MEN) = WNP1(1,MEN)*Kn**1.5d0/2.44/DW(MEN)*MR
    end if
  
    !Part 10:
    St(1,ME)=0d0
    St(1,ME)=0d0
          
  end do

!*******************************************************************************************
    End
!###########################################################################################


    
    

    
    