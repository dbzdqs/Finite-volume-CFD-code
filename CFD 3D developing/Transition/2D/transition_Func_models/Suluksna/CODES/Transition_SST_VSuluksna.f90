!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// Transition V Suluksna                                                                //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Transition_SST_VSuluksna(Dim,Ncyc,INW,X,Y,NX,NY,NC,NF,NF1,NF2,NFW1,NFW2,NFI1,&
                                     NFI2,NFO1,NFO2,NFS1,NFS2,NFF1,NFF2,NP,IDS,XC,YC,DW,DT,A,&
									 Minf,Rinf,MR,R0,NRKS,Mu0,Wb,Wbt,Wnp1,Wnt,Wntp1,Mu,Mut,&
									  NN,Mutinf,TUinf)
 Implicit None
!********************************************************************************************** 
 Integer::Dim,I,J,Ii,Jj,NC,NF,NFW1,NFW2,NF1,NF2,NS,Ncyc,NP,NRKS,ME,NFI1,NFI2,NFO1,NFO2,NFS1,&
          NFS2,NFF1,NFF2  ,NN
 Real(8)::ALF,Minf,Mu0,U,V,RKco,Rinf,K,Omega,R0,MR,Mutinf,TUinf
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:Dim)::X,Y,NX,NY,A,Mu,Mut,Mutq,XC,YC,DW,DT,DKX,DKY,DOmegX,DOmegY,DGX,DGY,&
                           DRX,DRY
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:4,1:Dim)::Wntp1,Wbt,Wnt,Cont,Dift,Dist,Sout,Wct
 Real(8),Dimension(1:Dim)::DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,DDGX,DDGY,DDRX,DDRY,&
                           CDKw,F11,F22,sigmak,sigmaw,beta,gama,Vor,Error
 Real(8)::sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,beta2,gama1,gama2,betastar,Rho
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
    Call TraSST_BC_V2(Dim,NX,NY,NF,NFS1,NFS2,NFO1,NFO2,NFW1,NFW2,NFI1,NFI2,NFF1,NFF2,DW,Wb,&
	                  Wbt,Wnp1,Wntp1,IDS,Minf,Rinf,MR,Mu,Mutinf,TUinf)

  !************************************************************************

 !Part 3:
    Do I=1,NC
       Wnt(1,I) =Wntp1(1,I)
       Wnt(2,I) =Wntp1(2,I)
       Wnt(3,I) =Wntp1(3,I)
       Wnt(4,I) =Wntp1(4,I)
       Mutq(I)  =Mut(I)
    End Do

      !Part 4:
      Do NS=1,1 !NRKS
          
      !Part 5:   
	   RKco=0.1 !/(NRKS-NS+1)  0.4
       
      !Part 6:
       Call TraSST_GradCell(A,Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DW,MR,Wnp1,Wntp1,Wb,Wbt,&
                           DDGX,DDGY,DDRX,DDRY,Mu,Mut)
       Call KwSST_GradCell(A,Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DW,MR,Wnp1,Wntp1,Wb,Wbt,&
	                       DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,DDOmegX,DDOmegY,&
                            Mu,Mut)
    
      !Part 7:
       Call TraSST_GradFace(Dim,NC,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,DGX,DGY,DRX,DRY,&
	                        DW,MR,beta1,Wnp1,Wntp1,Wb,Wbt,Mu,Mut)
       Call KwSST_GradFace(Dim,NC,NFW1,NFW2,NF,NF1,NF2,NP,IDS,X,Y,XC,YC,DKX,DKY,DOmegX,&
	                       DOmegY,DW,MR,beta1,Wnp1,Wntp1,Wb,Wbt,Mu,Mut)
       
       !part 8:
       Call KwSST_Funcs_V1(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,DDKX,&
	                       DDKY,DDOmegX,DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,beta1,&
						   beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,beta,gama,Mu,Mut)
                  
       
       Call TraSST_Funcs_VSuluksna(A,Dim,NC,IDS,NX,NY,DW,MR,Wnp1,Wntp1,DDUX,DDUY,DDVX,DDVY,&
	                               DDKX,DDKY,DDOmegX,DDOmegY,sigmak1,sigmak2,sigmaw1,sigmaw2,&
								   beta1,beta2,gama1,gama2,betastar,F11,F22,Sigmak,Sigmaw,beta,&
								   gama,Mu,Mut,Flength,Strain,Fonset,Fturb,T,Reteq,Ftt,Geff)

      !Part 9:
       Call KwSST_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,Wntp1,Wb,Wbt,Cont)
	   Call TraSST_Con(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,Wnp1,Wntp1,Wb,Wbt,Cont)
       
       
       
      !Part 10:
	   Call KwSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX,DKY,DOmegX,DOmegY,MR,Wnp1,&
	                  Wntp1,Wb,Wbt,Mu,Mut,Sigmak,Sigmaw,beta1,Dift)
       Call TraSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DGX,DGY,DRX,DRY,MR,Wnp1,Wntp1,&
	                   Wb,Wbt,Mu,Mut,Dift)
	               
      !Part 11:
       Call TraSST_Source(A,Dim,IDS,NC,MR,Wnp1,Wntp1,Mu,Mut,DDUX,DDUY,DDVX,DDVY,DDKX,DDKY,&
	                      DDOmegX,DDOmegY,betastar,sigmaw2,beta,gama,F11,Flength,Strain,&
						  Fonset,Fturb,T,Reteq,Ftt,Geff,Sout)
      
      !Part 12:
      Do J=1,NC
          
          Wntp1(1,J) = Wnt(1,J)-RKco*DT(J)*( Cont(1,J)/A(J)+Dift(1,J)/A(J) + Sout(1,J)  ) 
          Wntp1(2,J) = Wnt(2,J)-RKco*DT(J)*( Cont(2,J)/A(J)+Dift(2,J)/A(J) + Sout(2,J)  ) 
          Wntp1(3,J) = Wnt(3,J)-RKco*DT(J)*( Cont(3,J)/A(J)+Dift(3,J)/A(J) + Sout(3,J)  ) 
          Wntp1(4,J) = Wnt(4,J)-RKco*DT(J)*( Cont(4,J)/A(J)+Dift(4,J)/A(J) + Sout(4,J)  )

                    
      !Part 13: 
          if (Wntp1(1,J)<0.0 ) then
           !   print*, "K" 
              Wntp1(1,J)=Wnt(1,J)
          end if

          if (Wntp1(2,J)<0.0 ) then
            !  print*, "Omega"
             Wntp1(2,J)=Wnt(2,J)
          end if
          
         if (Wntp1(3,J)<0.0 .OR. Wntp1(3,J)>1.0 ) then
            !  print*, "Gama" 
              Wntp1(3,J)=Wnt(3,J)
         end if
         
          if (Wntp1(4,J)<0.0 ) then
            !  print*, "Re Theta" 
              Wntp1(4,J)=Wnt(4,J)
          end if           



      !Part 14:
          K = Wntp1(1,J)/Wnp1(1,J)
          Omega = Wntp1(2,J)/Wnp1(1,J)
          Rho = Wnp1(1,J)
                    
      !Part 15:
          Vor(J)= Dsqrt ( (DDUY(J)-DDVX(J))*(DDUY(J)-DDVX(J)) )
          
      !Part 16:     
          Mut(J) = (0.31*Rinf/Minf)*Rho*K / max(0.31*Omega,Vor(J)*F22(J)) 
              
      End Do !J       

      END DO !NS
                            
    

!****************************************************************************************      
      !Part 17:
      Do I=1,NC
          Error(I) = Dabs(Mut(I)-Mutq(I))
      End Do
      
      if ( Mod(Ncyc,50)==0 ) then
      print*, "Error=" , maxval(Error(1:NC))
      end if
    
!****************************************************************************************
 End
!##############################################################################################
