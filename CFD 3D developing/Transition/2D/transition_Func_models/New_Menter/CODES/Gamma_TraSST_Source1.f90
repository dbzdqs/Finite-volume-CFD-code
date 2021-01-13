!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// calculation of Source term for transition model                                      //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Gamma_TraSST_Source1(A,Dim,IDS,NC,MR,Wnp1,Wntp1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,&
                                 DDKY,DDOmegX,DDOmegY,betastar,sigmaw2,beta,gama,F11,Strain,&
								 Fonset,Fturb,T,Reteq,Ftt,Geff,Rev,Sout)
!**********************************************************************************************
 Integer::Dim,I,NC,J,ME
 Real(8)::MR,betastar,sigmaw2,K,Omega,Rho,Gamma,Rett,Dk
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,toxx,toxy,toyy,&
                           prow,Prok,beta,gama,F11,A,Vor,Pbark,Dbark,Pk,Pgama,Egama,Plimk,&
						   Flimon,Rev
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Sout
 Real(8),Dimension(1:Dim)::Strain,Fonset,Fturb,T,Reteq,Ftt,Geff,Modi
 Real(8)::ce1,ce2,ca1,ca2,ctt,Flength,sigmaG,ck,Csep,Relimtc
!***********************************************************************************************
 !Part1! Constants
 Flength=100.0
 ce2 = 50.0
 ca2 = 0.06
 sigmaG=1.0
 ck=1.0
 Csep=1.0
 Relimtc=1100.0
 
 
 !Part 2:
 Do I=1,NC
     Sout(1,I)=0.0
     Sout(2,I)=0.0
     Sout(3,I)=0.0
 End Do
  
 !Part 3:
 Do I=1,NC
     
    !Part 4: 
    Rho=Wnp1(1,I)
    K=Wntp1(1,I)/Wnp1(1,I)
    Omega=Wntp1(2,I)/Wnp1(1,I)
    Gamma = Wntp1(3,I)/Wnp1(1,I)
    
    !Part 5:
    toxx(I) = (MR)*Mut(I)*( (4.0/3.0)*DDUX(I)-(2.0/3.0)*DDVY(I)  ) -(2.0/3.0)*Rho*K
    toxy(I) = (MR)*Mut(I)*( DDUY(I)+DDVX(I)  )
    toyy(I) = (MR)*Mut(I)*( (4.0/3.0)*DDVY(I)-(2.0/3.0)*DDUX(I)  ) -(2.0/3.0)*Rho*K
    Vor(I) =  Dsqrt ( (DDUY(I)-DDVX(I))*(DDUY(I)-DDVX(I)) )
    Strain(I) = Dsqrt(2.0*DDUX(I)*DDUX(I) + (DDUY(I)+DDVX(I))*(DDUY(I)+DDVX(I)) +&
	            2.0*(DDVY(I)*DDVY(I)))   
    
    !Part 6:
    Prow(I) = toxx(I)*DDUX(I) + toxy(I)*(DDUY(I)+DDVX(I)) + toyy(I)*DDVY(I)    
    Prow(I) = abs(Prow(I))

    !Part 7:
    Prok(I) = min (Prow(I),20.0*betastar*Rho*Omega*K)
    Modi(I) = max(Gamma,0.1)

    Flimon(I)=min(max(Rev(I)/(2.2*Relimtc)-1.0,0.0),3.0)
    Plimk(I)=5.0*ck*max(Gamma-0.2,0.0)*(1.0-Gamma)*Flimon(I)*max(3.0*Csep*Mu(I)*&
	         MR-Mut(I)*MR,0.0)*Strain(I)*Vor(I)
    
    Pgama(I)=Flength*Rho*Strain(I)*Gamma*(1.0-Gamma)*Fonset(I)
    Egama(I)=ca2*Rho*Vor(I)*Gamma*Fturb(I)*(ce2*Gamma-1.0)
    
    
    !Part 8:
    Sout(1,I) = ( -Gamma*Prok(I) + (betastar*Rho*Omega*K)*Modi(I) - Plimk(I) )   
    Sout(2,I) = ( -(1.0/MR)*(gama(I)*Rho/Mut(I))*Prow(I) + beta(I)*Rho*Omega*Omega &
    - 2.0*(1.0-F11(I))*Rho*sigmaw2*(DDKX(I)*DDomegX(I)+DDKY(I)*DDomegY(I))/Omega  )
    
    Sout(3,I) = ( -Pgama(I) + Egama(I) )
    

 END DO
 
!********************************************************************************************** 
End
!##############################################################################################

