!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!// calculation of diffiusion term for transition model                                  //!
!// Date: October, 12, 2017                                                              //!
!// Developed by: M.A.Zoljanahi, Iran, Tehran, OpenFlows@chmail.ir                       //!
!// Doc ID: MC2F008F1                                                                    //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Subroutine Gamma_TraSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DGX,DGY,MR,Wnp1,Wntp1,Wb,&
                             Wbt,Mu,Mut,Dift)
 Implicit None
!**********************************************************************************************
 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE
 Real(8)::K,Omeg,DX,DY,AREA,DTX,DTY,MR,Mum,Mutm,Sigmakm,Sigmawm,beta1
 Real(8),Dimension(1:4,1:Dim)::Wnp1
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:Dim)::Mu,Mut,Sigmak,Sigmaw,NX,NY,DKX,DKY,DOmegX,DOmegY,DGX,DGY
 Real(8),Dimension(1:3,1:Dim)::Wntp1,Wbt,Dift
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NEC
 Real(8)::sigmaf,sigmatt
!***************************************************************************************************
 !Constants:
 sigmaf =  1.0
 sigmatt = 2.0
 
 
 !Part 1:
DO I=1,NC
    Dift(3,I) = 0.0
END DO

 
    !Part 2:
 DO I=NF1+1,NF2
    
    !Part 3: 
    ME = IDS(1,I)
    NE = IDS(2,I)

    !Part 4:
    Mum = 0.5*(Mu(ME)+Mu(NE))
    Mutm = 0.5*(Mut(ME)+Mut(NE))
   
    
   !Part 5:
    Dift(3,ME) = Dift(3,ME) + ( (Mum+Mutm*Sigmaf)*DGX(I)*NX(I) + (Mum+Mutm*Sigmaf)*DGY(I)*NY(I) )
    
    Dift(3,NE) = Dift(3,NE) - ( (Mum+Mutm*Sigmaf)*DGX(I)*NX(I) + (Mum+Mutm*Sigmaf)*DGY(I)*NY(I) )

 END DO

    !Part 6:
 DO I=NFW1+1,NFW2
  
   !Part 7:
    ME = IDS(1,I)
    
   !Part 8:
    Mum =  Mu(ME)
    Mutm = 0.0
      
   !Part 9:
    Dift(3,ME) = Dift(3,ME) + (Mum+Mutm*Sigmaf)*DGX(I)*NX(I) + (Mum+Mutm*Sigmaf)*DGY(I)*NY(I)
    
 END DO
 
   !Part 10:
   DO I=NFW2+1,NF
  
   !Part 11:
    ME = IDS(1,I)
  
   !Part 12:
    Mum =  Mu(ME)
    Mutm = Mut(ME)
      
   !Part 13:
    Dift(3,ME) = Dift(3,ME) + (Mum+Mutm*Sigmaf)*DGX(I)*NX(I) + (Mum+Mutm*Sigmaf)*DGY(I)*NY(I)
     
  END DO

 
!Part 14:
 Do I=1,NC
    
    Dift(3,I) = -(MR)*Dift(3,I)

 End do
 

!**********************************************************************************************
 End
!##############################################################################################

