!!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Main Transition 2D subroutine                                           //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!**********************************************************************************************
 Subroutine Gamma_Transition_SST(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,NFI2,NFO1,&
                                 NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,Minf,Rinf,MR,R0,&
                                 NRKS,Mu0,Wb,Wbt,Wnp1,Wnt,Wntp1,Mu,Mut,TUinf,Mutinf)
 Implicit None
!********************************************************************************************** 
 Integer::Dim,I,J,Ii,Jj,NC,NF,NFW1,NFW2,NF1,NF2,NS,Ncyc,NP,NRKS,ME,NFI1,NFI2,NFO1,NFO2,NFS1,NFS2,&
          NFF1,NFF2  ,NN
 Real(8)::ALF,Minf,Mu0,U,V,RKco,Rinf,K,Omega,R0,MR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,A,Mu,Mut,Mutq,XC,YC,DW,DT,DKX,DKY,DOmegX,DOmegY,DGX,DGY,DRX,DRY
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Wbt,Wnt,Cont,Dift,Dist,Sout,Wct
 
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,DDGX,DDGY,DDRX,DDRY,CDKw,&
                           F11,F22,sigmak,sigmaw,beta,gama,Vor,Error,Rev
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,Rho,Mutinf,TUinf
 
 Real(8),Dimension(1:Dim)::Flength,Strain,Fonset,Fturb,T,Reteq,Ftt,Geff

!**********************************************************************************************
    !Part 1:
    
    sigmak1=0.85
    sigmak2=1.0
    sigmaw1=0.5
    sigmaw2=0.856
    beta1=0.075
    beta2=0.0828
    gama1=5.0/9.0
    gama2=0.44
    betastar=0.09
        
    !Part 2:
    Call Gamma_TraSST_BC(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,Wbt,&
	                     Wnp1,Wntp1,IDS,Minf,Rinf,MR,Mu,TUinf,Mutinf)

  !************************************************************************

 !Part 3:
    Do I=1,NC
       Wnt(1,I) =Wntp1(1,I)
       Wnt(2,I) =Wntp1(2,I)
       Wnt(3,I) =Wntp1(3,I)
       Mutq(I)  =Mut(I)
    End Do

      !Part 4:
      Do NS=1,NRKS
              
      !Part 5:   
	   RKco=1.0/(NRKS-NS+1.0)

      !Part 6:
       Call Gamma_KwSST_GradCell(A,Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DW,MR,Wnp1,Wntp1,Wb,Wbt,DDUX,&
	                             DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,Mu,Mut)
    
      !Part 7:
       Call Gamma_TraSST_GradFace(Dim,NC,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,DGX,DGY,DW,MR,beta1,Wnp1,&
	                              Wntp1,Wb,Wbt,Mu,Mut)
       Call Gamma_KwSST_GradFace(Dim,NC,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,DKX,DKY,DOmegX,DOmegY,DW,MR,&
	                             beta1,Wnp1,Wntp1,Wb,Wbt,Mu,Mut)
       
       !part 8:
       Call Gamma_KwSST_Funcs(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,&
	                          DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,F11,&
							  F22,Sigmak,Sigmaw,beta,gama,Mu,Mut)
                  
       
       Call Gamma_TraSST_Funcs(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,&
	                           DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,&
							   F11,F22,Sigmak,Sigmaw,beta,gama,Mu,Mut,Flength,Strain,Fonset,Fturb,T,Reteq,&
							   Ftt,Rev)

      !Part 9:
       Call Gamma_KwSST_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,Wntp1,Wb,Wbt,Cont)
	   Call Gamma_TraSST_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,Wntp1,Wb,Wbt,Cont)
       
       
       
      !Part 10:
	   Call Gamma_KwSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX,DKY,DOmegX,DOmegY,MR,Wnp1,Wntp1,Wb,Wbt,&
	                        Mu,Mut,Sigmak,Sigmaw,beta1,Dift)
       Call Gamma_TraSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DGX,DGY,MR,Wnp1,Wntp1,Wb,Wbt,Mu,Mut,Dift)
	               
      !Part 11:
       Call Gamma_TraSST_Source1(A,Dim,IDS,NC,MR,Wnp1,Wntp1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,&
	                             DDOmegY,betastar,sigmaw2,beta,gama,F11,Strain,Fonset,Fturb,T,Reteq,Ftt,Geff,&
								 Rev,Sout)
       
      !Part 12:
      Do J=1,NC
          
          Wntp1(1,J) = Wnt(1,J)-RKco*DT(J)*( Cont(1,J)/A(J)+Dift(1,J)/A(J) + Sout(1,J)  ) 
          Wntp1(2,J) = Wnt(2,J)-RKco*DT(J)*( Cont(2,J)/A(J)+Dift(2,J)/A(J) + Sout(2,J)  ) 
          Wntp1(3,J) = Wnt(3,J)-RKco*DT(J)*( Cont(3,J)/A(J)+Dift(3,J)/A(J) + Sout(3,J)  ) 

                    
      !Part 13: 
          if (Wntp1(1,J)<0.0 ) then
              Wntp1(1,J)=Wnt(1,J)
          end if

          if (Wntp1(2,J)<0.0 ) then
             Wntp1(2,J)=Wnt(2,J)
          end if
          
         if (Wntp1(3,J)<0.0 .OR. Wntp1(3,J)>1.0 ) then
              Wntp1(3,J)=Wnt(3,J)
         end if
                



      !Part 14:
          K = Wntp1(1,J)/Wnp1(1,J)
          Omega = Wntp1(2,J)/Wnp1(1,J)
          Rho = Wnp1(1,J)
                    
      !Part 15:
          Vor(J)= Dsqrt ( (DDUY(J)-DDVX(J))*(DDUY(J)-DDVX(J)) )
          Strain(I) = Dsqrt(2.0*DDUX(I)*DDUX(I) + (DDUY(I)+DDVX(I))*(DDUY(I)+DDVX(I)) +&
		                    2.0*(DDVY(I)*DDVY(I)))
   
      !Part 16:     
          Mut(J) = (0.31)*Rho*K / max(0.31*Omega, Strain(I)*F22(J)) / MR
              
      End Do       

      END DO 
                                 
      !Part 17:
      Do I=1,NC
          Error(I) = Dabs(Mut(I)-Mutq(I))
      End Do
      
      if ( Mod(Ncyc,10)==0 ) then
      print*, "Error=" , maxval(Error(1:NC))
      end if
 
!****************************************************************************************
 End
!##############################################################################################
