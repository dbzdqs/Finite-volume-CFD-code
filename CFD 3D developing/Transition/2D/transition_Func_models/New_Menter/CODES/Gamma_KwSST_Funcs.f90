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
!*******************************************************************************************
 Subroutine Gamma_KwSST_Funcs(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,&
                              DDKY,DDOmegX,DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,&
							  beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,beta,gama,&
							  Mu,Mut)
                  
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
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,k,Omega,Rho,u,v,&
          Rett,Gamma
 
 Real(8),Dimension(1:Dim)::Retc,Strain,Rev,Fonset1,Fonset2,Fonset3,Fonset,RT,Fturb,Rew,&
                           Fsublayer,Flength1,Flength,UU,T,Delta,Fwake,Ftt,TU,dU,tt,Lt,AA,&
						   FLt,FLtp,Reteq,Reteqp,fx,fpx,ttn1,Freattach,Gsep,Geff,&
                           Ry,F33
 Real(8)::ce2
!***************************************************************************************************  
! Functions

!Part 1:
Do I=1,NC
        
    k=Wntp1(1,I)/Wnp1(1,I)
    Omega=Wntp1(2,I)/Wnp1(1,I)
    Rho=Wnp1(1,I)
        
    
    CDkw(I)=max( 2.0*Rho*sigmaw2*(DDKX(I)*DDomegX(I)+DDKY(I)*DDomegY(I))/Omega,10.0**-10.0)
    
    arg1(I)=min (max (dsqrt(abs(K))/betastar/Omega/DW(I),500.0*Mu(I)*MR/DW(I)/DW(I)/Omega/Rho ), &
    4.0*Rho*sigmaw2*K/DW(I)/DW(I)/CDkw(I) )
    F11(I)=tanh(arg1(I)*arg1(I)*arg1(I)*arg1(I))
        
    arg2(I)=max (2.0*dsqrt(abs(K))/betastar/DW(I)/Omega ,500.0*Mu(I)*MR/DW(I)/DW(I)/Omega/Rho )
    F22(I)=tanh(arg2(I)*arg2(I))

        
  
    Ry(I) = Rho*DW(I)*dsqrt(abs(K))/(Mu(I)*MR)
    F33(I) = exp(-(Ry(I)/120.0)**8.0)
    F11(I) = max(F11(I) , F33(I))
    

    
    sigmak(I)=F11(I)*sigmak1 + (1.0-F11(I))*sigmak2
    sigmaw(I)=F11(I)*sigmaw1 + (1.0-F11(I))*sigmaw2
    beta(I)=F11(I)*beta1 + (1.0-F11(I))*beta2
    gama(I)=F11(I)*gama1 + (1.0-F11(I))*gama2
END Do
    
    

    
    End
!##############################################################################################

