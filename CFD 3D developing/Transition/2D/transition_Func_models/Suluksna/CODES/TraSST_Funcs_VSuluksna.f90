!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// KW_SST _Funcs Subroutine                                                             //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine TraSST_Funcs_VSuluksna(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,&
                                   DDKX,DDKY,DDOmegX,DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,&
								   beta1,beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,&
								   beta,gama,Mu,Mut,Flength,Strain,Fonset,Fturb,T,Reteq,Ftt,Geff)
 Implicit None
!**********************************************************************************************
 Integer::Dim,I,NC,IJ,ME
 Real(8)::AREA,MR
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut
 Real(8),Dimension(1:4,1:Dim)::Wntp1
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,NX,NY,CDKw,F11,F22,&
                           arg1,arg2,sigmak,sigmaw,beta,gama,Vor
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,k,Omega,Rho,u,v,&
          Rett,Gamma
 Real(8),Dimension(1:Dim)::Retc,Strain,Rev,Fonset1,Fonset2,Fonset3,Fonset,RT,Fturb,Rew,&
                           Fsublayer,Flength1,Flength,UU,T,Delta,Fwake,Ftt,TU,dU,tt,Lt,AA,&
						   FLt,FLtp,Reteq,Reteqp,fx,fpx,ttn1,Freattach,Gsep,Geff,Ry,F33,KB,FL,&
						   FKB,LP,FLP
 Real(8)::ce2
!***************************************************************************************************  
! Functions
      ce2 = 50.0

  !Part 1:
    
    
    !Part 2:
Do I=1,NC
     
    Rho=Wnp1(1,I)
    u=Wnp1(2,I)/Wnp1(1,I)
    v=Wnp1(3,I)/Wnp1(1,I)
    K=Wntp1(1,I)/Wnp1(1,I)
    Omega=Wntp1(2,I)/Wnp1(1,I)
    Gamma = Wntp1(3,I)/Wnp1(1,I)
    Rett = Wntp1(4,I)/Wnp1(1,I)
    
 
    !part 3:
    Retc(I)= min(max(-(0.025*Rett)**2 + 1.47*Rett - 120.0 , 125.0 ) , Rett )
   
    ! part 4:
    Strain(I) =  Dsqrt(2.0*DDUX(I)*DDUX(I) + (DDUY(I)+DDVX(I))*(DDUY(I)+DDVX(I)) + &
	                   2.0*(DDVY(I)*DDVY(I)))   
    Rev(I) = Rho*Strain(I)*DW(I)*DW(I)/ (Mu(I)*MR)    
    Fonset1(I) = Rev(I) / (2.193*Retc(I))   
    Fonset2(I) = min(max(Fonset1(I),Fonset1(I)**4.0),2.0)

    RT(I) = Rho*k / (MR*Mu(I)*Omega)
    Fonset3(I) = max( 1.0 - (RT(I)/2.5)**3.0 , 0.0)
    
    Fonset(I) = max(Fonset2(I)-Fonset3(I) , 0.0)

    Fturb(I) = exp(-(RT(I)/4.0)**4.0)
        
    Rew(I) = Rho*Omega*Dw(I)*Dw(I) / (MR*Mu(I))
    
    Flength(I) = min(0.1*exp(-0.022*Rett + 12.0 ) + 0.45 , 300.0 )
    

    
    !part 5:
    UU(I) = Dsqrt(u*u + v*v)
    T(I) = MR*500.0*Mu(I) / (Rho*UU(I)*UU(I))
    Vor(I) =  Dsqrt ( (DDUY(I)-DDVX(I))*(DDUY(I)-DDVX(I)) )
    Delta(I) = 375.0*MR*Vor(I)*Mu(I)*Rett*DW(I)/(Rho*UU(I)*UU(I))
    Fwake(I) = exp(-(Rew(I)/100000.0)**2.0)
    Ftt(I) = min(max(Fwake(I)*exp(-(DW(I)/Delta(I))**4.0) , 1.0 - &
	          ((ce2*gamma-1.0)/(ce2-1.0))**2.0) , 1.0)
    
       
    TU(I) = 100.0*Dsqrt(2.0*abs(k)/3.0)/UU(I)
    
    dU(I) = ( ((u*u)/(UU(I)*UU(I))) * DDUX(I) + ((u*v)/(UU(I)*UU(I))) * DDUY(I) + &
	        ((u*v)/(UU(I)*UU(I))) * DDVX(I) + ((v*v)/(UU(I)*UU(I))) * DDVY(I) )
    KB(I)=Mu(I)*MR*dU(I)/(Rho*UU(I)*UU(I))
    
    !part 6:
    ! Initial Guess
    Reteq(I)=803.73*(Tu(I)+0.6067)**(-1.027)
    
    if (Reteq(I)<20.0) then
        Reteq(I)=20.0
    end if

    tt(I) =  Reteq(I)*Mu(I)*MR/(Rho*UU(I))
    !!!!!!!!!
    
    Do ij=1,3
    AA(I) = Rho*dU(I)/(Mu(I)*MR)
    Lt(I) = AA(I)*tt(I)*tt(I)
    
    if (Lt(I)>=0.1) then
        Lt(I)=0.1
        AA(I)=0.1/tt(I)/tt(I)
    end if
    if (Lt(I)<=-0.1) then
        Lt(I)=-0.1
        AA(I)=-0.1/tt(I)/tt(I)
    end if
    
    FL(I) = 10.32*Lt(I) + 89.47*Lt(I)*Lt(I) + 265.51*Lt(I)*Lt(I)*Lt(I)
    FKB(I) = 0.0962*(KB(I)*10**6) + 0.148*(KB(I)*10**6)*(KB(I)*10**6) + &
	         0.0141*(KB(I)*10**6)*(KB(I)*10**6)*(KB(I)*10**6)
    LP(I) = 2*AA(I)*tt(I)
    FLP(I) = 10.32*LP(I) + 89.47*2*Lt(I)*LP(I) + 265.51*3*Lt(I)*Lt(I)*LP(I)
    
    
    if (Lt(I)<=0.0) then
        FLt(I) = 1 + FL(I)*exp(-Tu(I)/3.0)
        FLtp(I) = exp(-Tu(I)/3.0)*FLP(I)
    else
        FLt(I) = 1 + FKB(I)*(1-exp(-2.0*Tu(I)/3.0)) + 0.556*(1-exp(-23.9*Lt(I)))*&
		         exp(-Tu(I)/3.0)
        FLtp(I) = 0.556*exp(-Tu(I)/3.0)*(-1.0*exp(-23.9*Lt(I))*(-23.9)*LP(I))
    end if
    
        Reteq(I)= 803.73*((Tu(I)+0.6067)**(-1.027))*FLt(I)
        Reteqp(I) = 803.73*((Tu(I)+0.6067)**(-1.027))*FLtp(I)
        
    
    fx(I) = Reteq(I) - Rho*UU(I)*tt(I)/(Mu(I)*MR)
    fpx(I) = Reteqp(I) - Rho*UU(I)/(Mu(I)*MR)
    ttn1(I) = tt(I) - fx(I)/fpx(I)
    tt(I)=ttn1(I)
    end do

      
    !part 7:
    if (Reteq(I)<20.0) then
        Reteq(I)=20.0
    end if

   
    !part 8:
    Freattach(I) = exp(-(RT(I)/20.0)**4.0)
    Gsep(I) = min(2.0*Freattach(I)*max(0.0 , (Rev(I)/3.235/Retc(I))-1.0) , 2.0) * Ftt(I)
    Geff(I) = max(Gamma,Gsep(I))
    

END Do

!*********************************************************************************************    
    End
!##############################################################################################

