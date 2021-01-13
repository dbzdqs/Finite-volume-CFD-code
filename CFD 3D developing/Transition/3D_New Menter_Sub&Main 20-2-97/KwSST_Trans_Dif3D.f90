!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Diffusion Terms of Transition Model                       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M. Namvar,M.A.Zoljanahi Iran, Tehran, OpenFlows@chmail.ir              //!
!// Doc ID: MC2F058F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine KwSST_Trans_Dif3D(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,&
                            DOmegX_F,DOmegY_F,DOmegZ_F,DGamaX_F,DGamaY_F,DGamaZ_F,&
                            MR,Sigk,Sigw,Sigg,WNP1,WTNP1,Mu,Mut,Dift)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,&
                DOmegZ_F,DGamaX_F,DGamaY_F,DGamaZ_F,MR,WNP1,WTNP1,Mu,Mut,Sigk,Sigw,Sigg
                
 Intent(Out  )::Dift

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE
 Real(8)::Mu_k,Mu_w,Mu_g,MR,F1,F2,F3,sigmakm,sigmawm,Sigmagm,Mum,Mutm
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:3,1:Dim)::WTNP1,Dift
 Real(8),Dimension(1:Dim)::Mu,Mut,Sigk,Sigw,Sigg,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DOmegX_F,DOmegY_F,&
                           DOmegZ_F,DGamaX_F,DGamaY_F,DGamaZ_F
!*******************************************************************************************
!Part 1:
 Do I=1,NC
    Dift(1,I) = 0.0
    Dift(2,I) = 0.0
    Dift(3,I) = 0.0
 End do

!Part 2:
 Do I=NF1+1,NF2
    
   !Part 3: 
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 4:
    Mum     = 0.5*(Mu(ME)+Mu(NE))
    Mutm    = 0.5*(Mut(ME)+Mut(NE))
    sigmakm = 0.5*(sigk(ME)+sigk(NE))
	sigmawm = 0.5*(sigw(ME)+sigw(NE))
    sigmagm = 0.5*(sigg(ME)+sigg(NE))
    
   !Part 5:
    Mu_k = Mum+Mutm*Sigmakm
	Mu_w = Mum+Mutm*Sigmawm
	Mu_g = Mum+Mutm*Sigmagm
    
	F1 = Mu_k * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) + DKZ_F(I)   *NZ(I) ) 
	F2 = Mu_w * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) + DOmegZ_F(I)*NZ(I) ) 
	F3 = Mu_g * ( DGamaX_F(I)*NX(I) + DGamaY_F(I)*NY(I) + DGamaZ_F(I)*NZ(I) )

   !Part 6:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    Dift(3,ME) = Dift(3,ME) + F3
    
    Dift(1,NE) = Dift(1,NE) - F1
    Dift(2,NE) = Dift(2,NE) - F2
    Dift(3,NE) = Dift(3,NE) - F3
    
 End do

!Part 7:
 Do I=NFW1+1,NFW2
  
   !Part 8:
    ME = IDS(1,I)
    
   !Part 9:    
	F1 = Mu(ME) * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) + DKZ_F(I)   *NZ(I) ) 
	F2 = Mu(ME) * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) + DOmegZ_F(I)*NZ(I) ) 
	F3 = Mu(ME) * ( DGamaX_F(I)*NX(I) + DGamaY_F(I)*NY(I) + DGamaZ_F(I)*NZ(I) )

   !Part 10:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    Dift(3,ME) = Dift(3,ME) + F3
    
 End do
 
!Part 11:
 Do I=NFW2+1,NF
  
   !Part 12:
    ME = IDS(1,I)
  
   !Part 13:
    Mu_k = Mu(ME) + Mut(ME)*Sigk(ME)
	Mu_w = Mu(ME) + Mut(ME)*Sigw(ME)
	Mu_g = Mu(ME) + Mut(ME)*Sigg(ME)
    
	F1 = Mu_k * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) + DKZ_F(I)   *NZ(I) ) 
	F2 = Mu_w * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) + DOmegZ_F(I)*NZ(I) )
	F3 = Mu_g * ( DGamaX_F(I)*NX(I) + DGamaY_F(I)*NY(I) + DGamaZ_F(I)*NZ(I) ) 

   !Part 14:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    Dift(3,ME) = Dift(3,ME) + F3
 
 End do

!Part 15:
 Do I=1,NC
    Dift(1,I) = -MR*Dift(1,I)
	Dift(2,I) = -MR*Dift(2,I)
	Dift(3,I) = -MR*Dift(3,I)
 End do
!*******************************************************************************************
 End
!###########################################################################################

