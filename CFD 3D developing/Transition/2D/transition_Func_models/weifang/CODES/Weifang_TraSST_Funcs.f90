!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// Transition V Weifang                                                                 //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Weifang_TraSST_Funcs(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,&
                  sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,beta,gama,Mu,Mut,&
                  Flength,Strain,Fonset,Fturb,T,Reteq,Ftt,Geff)
 Implicit None
!**********************************************************************************************
 Integer::Dim,I,NC,IJ,ME
 Real(8)::AREA,MR
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut
 Real(8),Dimension(1:3,1:Dim)::Wntp1
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,NX,NY,CDKw,F11,F22,arg1,arg2,sigmak,sigmaw,beta,gama,Vor
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,k,Omega,Rho,u,v,Rett,Gamma
 
 Real(8),Dimension(1:Dim)::Retc,Strain,Rev,Fonset1,Fonset2,Fonset3,Fonset,RT,Fturb,Rew,Fsublayer,Flength1,Flength,&
                           UU,T,Delta,Fwake,Ftt,TU,dU,tt,Lt,AA,FLt,FLtp,Reteq,Reteqp,fx,fpx,ttn1&
                           ,Freattach,Gsep,Geff,&
                            Ry,F33,T1
 Real(8)::ce2
!***************************************************************************************************  
! Functions
   !Part 1:
   ce2 = 50.0

 

    
    
    !Part 2:
Do I=1,NC
     
    Rho=Wnp1(1,I)
    u=Wnp1(2,I)/Wnp1(1,I)
    v=Wnp1(3,I)/Wnp1(1,I)
    K=Wntp1(1,I)/Wnp1(1,I)
    Omega=Wntp1(2,I)/Wnp1(1,I)
    Gamma = Wntp1(3,I)/Wnp1(1,I)
    
    
    !-///////////////////////////////////////////////////////Retc(I)///////
    !part 3:
    Vor(I) =  Dsqrt ( (DDUY(I)-DDVX(I))*(DDUY(I)-DDVX(I)) )
    Strain(I) =  Dsqrt(2.0*DDUX(I)*DDUX(I) + (DDUY(I)+DDVX(I))*(DDUY(I)+DDVX(I)) + 2.0*(DDVY(I)*DDVY(I)))
    RT(I) = Rho*k / (MR*Mu(I)*Omega)
    T1(I) = RT(I)*Vor(I)/Omega
    Retc(I) = 900.0 - 810.0*min(T1(I)*(0.285*Tu(I)**0.526 + 0.35 ) / 4.0 , 1.0 )
   
    ! part 4:   
    Rev(I) = Rho*Strain(I)*DW(I)*DW(I)/ (Mu(I)*MR)    
    Fonset1(I) = Rev(I) / (2.193*Retc(I))   
    Fonset2(I) = min(max(Fonset1(I),Fonset1(I)**4.0),2.0)

    
    Fonset3(I) = max( 1.0 - (RT(I)/2.5)**3.0 , 0.0)
    
    Fonset(I) = max(Fonset2(I)-Fonset3(I) , 0.0)

    Fturb(I) = exp(-(RT(I)/4.0)**4.0)
        
    Rew(I) = Rho*Omega*Dw(I)*Dw(I) / (MR*Mu(I))
    
    Flength(I) = max(0.1,30.0*Dlog(Tu(I)) + 89.97 )
    !///////////////////////////////////////////////////////////////////////////////////////////
    
    !part 5:
    UU(I) = Dsqrt(u*u + v*v)
    T(I) = MR*500.0*Mu(I) / (Rho*UU(I)*UU(I))
    Delta(I) = 375.0*MR*Vor(I)*Mu(I)*Rett*DW(I)/(Rho*UU(I)*UU(I))
    Fwake(I) = exp(-(Rew(I)/100000.0)**2.0)
    Ftt(I) = min(max(Fwake(I)*exp(-(DW(I)/Delta(I))**4.0) , 1.0 - ((ce2*gamma-1.0)/(ce2-1.0))**2.0) , 1.0)
    
       
    TU(I) = 100.0*Dsqrt(2.0*abs(k)/3.0)/UU(I)
    if (TU(I)<0.027) then 
        TU(I)=0.027
    end if
    dU(I) = ( ((u*u)/(UU(I)*UU(I))) * DDUX(I) + ((u*v)/(UU(I)*UU(I))) * DDUY(I) + ((u*v)/(UU(I)*UU(I))) * DDVX(I) + ((v*v)/(UU(I)*UU(I))) * DDVY(I) )
    
    
    
   

    
    !///////////////////////////////////////
    !part 6:
    Freattach(I) = exp(-(RT(I)/20.0)**4.0)
    Gsep(I) = min(2.0*Freattach(I)*max(0.0 , (Rev(I)/3.235/Retc(I))-1.0) , 2.0) * Ftt(I)
    Geff(I) = max(Gamma,Gsep(I))
    

END Do

    
    End
!##############################################################################################

