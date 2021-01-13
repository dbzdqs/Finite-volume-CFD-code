!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: Calculate the Diffusion Terms of Turbulence Model                       //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: October, 12, 2014                                                              //!
!// Developed by: M. Namvar,M.H Saadat Iran, Tehran, OpenFlows@chmail.ir                 //!
!// Doc ID: MC2F048F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Ke_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,&
                        MR,Mu,Mut,Dift)
 Implicit None
!*******************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F,Sigk,Sige,MR,Mu,Mut
 Intent(Out  )::Dift

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE
 Real(8)::Mu_k,Mu_w,Sigk,Sige,MR,F1,F2
 Real(8),Dimension(1:Dim)::Mu,Mut,NX,NY,DKX_F,DKY_F,DEPSX_F,DEPSY_F
 Real(8),Dimension(1:2,1:Dim)::Dift
 Integer,Dimension(1:4,1:Dim)::IDS
!*******************************************************************************************
!Part 1:
 Do I=1,NC
    Dift(1,I) = 0.0
    Dift(2,I) = 0.0
 End do

!Part 2:
 Do I=NF1+1,NF2
    
   !Part 3: 
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 4:
    Mu_k = 0.5*( Mu(ME) + Mut(ME)*Sigk + Mu(NE) + Mut(NE)*Sigk )
	Mu_w = 0.5*( Mu(ME) + Mut(ME)*Sige + Mu(NE) + Mut(NE)*Sige )
    
	F1 = Mu_k * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) ) 
	F2 = Mu_w * ( DEPSX_F(I)*NX(I) + DEPSY_F(I)*NY(I) ) 

   !Part 5:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    
    Dift(1,NE) = Dift(1,NE) - F1
    Dift(2,NE) = Dift(2,NE) - F2
    
 End do

!Part 6:
 Do I=NFW1+1,NFW2
  
   !Part 7:
    ME = IDS(1,I)
    
   !Part 8:    
	F1 = Mu(ME) * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) ) 
	F2 = Mu(ME) * ( DEPSX_F(I)*NX(I) + DEPSY_F(I)*NY(I) ) 

   !Part 9:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    
 End do
 
!Part 10:
 Do I=NFW2+1,NF

    ME = IDS(1,I)
  
    Mu_k = Mu(ME) + Mut(ME)*Sigk
	Mu_w = Mu(ME) + Mut(ME)*Sige
    
	F1 = Mu_k * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) ) 
	F2 = Mu_w * ( DEPSX_F(I)*NX(I) + DEPSY_F(I)*NY(I) ) 

    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
 
 End do

!Part 13:
 Do I=1,NC
    Dift(1,I) = -MR*Dift(1,I)
	Dift(2,I) = -MR*Dift(2,I)
 End do
!*******************************************************************************************
 End
!###########################################################################################

