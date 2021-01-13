!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// calculation of Functions for transition model                                        //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Gamma_TraSST_Funcs(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,&
                               DDKY,DDOmegX,DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,&
							   beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,beta,gama,Mu,&
							   Mut,Flength,Strain,Fonset,Fturb,T,Reteq,Ftt,Rev)
 Implicit None
!**********************************************************************************************
 Integer::Dim,I,NC,IJ,ME
 Real(8)::AREA,MR 
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::A,DW,Mu,Mut
 Real(8),Dimension(1:3,1:Dim)::Wntp1
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,NX,NY,CDKw,F11,F22,&
                           arg1,arg2,sigmak,sigmaw,beta,gama,Vor
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,k,Omega,Rho,u,v,Gamma
 
 Real(8),Dimension(1:Dim)::Retc,Strain,Rev,Fonset1,Fonset2,Fonset3,Fonset,RT,Fturb,Rew,Fsublayer,&
                           Flength1,Flength,UU,T,Delta,Fwake,Ftt,dU,tt,Lt,AA,FLt,FLtp,Reteq,&
						   Reteqp,fx,fpx,ttn1,Freattach,Gsep,Ry,F33,TUL,FPG,Ltl
 Real(8)::ce2,CTU1,CTU2,CTU3,CPG1,ClimPG1,CPG2,CPG3,ClimPG2
!***************************************************************************************************  
! part 1:
 ce2 = 50.0
 
 CTU1=100.0
 CTU2=1000.0
 CTU3=1.0
 
 CPG1=14.68
 ClimPG1=1.5
 CPG2=-7.34
 CPG3=0.0
 ClimPG2=3.0
 
    
    !Part 2:
Do I=1,NC
     
    Rho=Wnp1(1,I)
    u=Wnp1(2,I)/Wnp1(1,I)
    v=Wnp1(3,I)/Wnp1(1,I)
    K=Wntp1(1,I)/Wnp1(1,I)
    Omega=Wntp1(2,I)/Wnp1(1,I)
    Gamma = Wntp1(3,I)/Wnp1(1,I)
    
    !part 3:
    UU(I) = Dsqrt(u*u + v*v)
    dU(I) = ( ((u*u)/(UU(I)*UU(I))) * DDUX(I) + ((u*v)/(UU(I)*UU(I))) * DDUY(I) + &
	        ((u*v)/(UU(I)*UU(I))) * DDVX(I) + ((v*v)/(UU(I)*UU(I))) * DDVY(I) )
    Ltl(I)=7.57*10**(-3)*dU(I)*Rho*Dw(I)*Dw(I)/(Mu(I)*MR)+0.0128
    
    Ltl(I)=min(max(Ltl(I),-1.0),1.0)
    
    if (Ltl(I)>=0.0) then
        FPG(I)=min(1+CPG1*Ltl(I),ClimPG1)
    end if
    if (Ltl(I)<0.0) then
        FPG(I)=min(1+CPG2*Ltl(I)+CPG3*min(Ltl(I)+0.0681,0.0),ClimPG2)
    end if
    
    FPG(I)=max(FPG(I),0.0)
    
    TUL(I)=min(100.0*Dsqrt(2.0*abs(k)/3.0)/(Omega*Dw(I)),100.0)
    Retc(I)=CTU1+CTU2*exp( -CTU3*TUL(I)*FPG(I) )

    ! part 4:
    Strain(I) =  Dsqrt(2.0*DDUX(I)*DDUX(I) + (DDUY(I)+DDVX(I))*(DDUY(I)+DDVX(I)) +&
	                   2.0*(DDVY(I)*DDVY(I)))   
    Rev(I) = Rho*Strain(I)*DW(I)*DW(I)/ (Mu(I)*MR)    
    Fonset1(I) = Rev(I) / (2.2*Retc(I))   
    Fonset2(I) = min(Fonset1(I),2.0)

    RT(I) = Rho*k / (MR*Mu(I)*Omega)
    Fonset3(I) = max( 1.0 - (RT(I)/3.5)**3.0 , 0.0)
    
    Fonset(I) = max(Fonset2(I)-Fonset3(I) , 0.0)

    Fturb(I) = exp(-(RT(I)/2.0)**4.0)
   

END Do

!**********************************************************************************************    
    End
!##############################################################################################

