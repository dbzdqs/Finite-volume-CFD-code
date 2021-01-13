!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: M. H. Saadat, Aerospace Eng., Amirkabir University of Technology         //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KwSST_Dif(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DOmegX_F,DOmegY_F,&
                      MR,Sigk,Sigw,WNP1,WTNP1,Mu,Mut,Dift)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,DKX_F,DKY_F,DOmegX_F,DOmegY_F,MR,WNP1,WTNP1,Mu,Mut,Sigk,Sigw
 Intent(Out  )::Dift

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE
 Real(8)::Mu_k,Mu_w,MR,F1,F2
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:2,1:Dim)::WTNP1,Dift
 Real(8),Dimension(1:Dim)::Mu,Mut,Sigk,Sigw,NX,NY,DKX_F,DKY_F,DOmegX_F,DOmegY_F
!*********************************************************************************************
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
    Mu_k = 0.5*( Mu(ME) + Mut(ME)*Sigk(ME) + Mu(NE) + Mut(NE)*Sigk(NE) )
	Mu_w = 0.5*( Mu(ME) + Mut(ME)*Sigw(ME) + Mu(NE) + Mut(NE)*Sigw(NE) )
    
	F1 = Mu_k * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) ) 
	F2 = Mu_w * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) ) 

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
	F1 = Mu(ME) * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) ) 
	F2 = Mu(ME) * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) ) 

   !Part 5:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    
 End do
 
!Part 10:
 Do I=NFW2+1,NF
  
   !Part 11:
    ME = IDS(1,I)
  
    Mu_k = Mu(ME) + Mut(ME)*Sigk(ME)
	Mu_w = Mu(ME) + Mut(ME)*Sigw(ME)
    
	F1 = Mu_k * ( DKX_F(I)   *NX(I) + DKY_F(I)   *NY(I) ) 
	F2 = Mu_w * ( DOmegX_F(I)*NX(I) + DOmegY_F(I)*NY(I) ) 

   !Part 5:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
 
 End do

!Part 14:
 Do I=1,NC
    Dift(1,I) = -MR*Dift(1,I)
	Dift(2,I) = -MR*Dift(2,I)
 End do
!*********************************************************************************************
 End
!###########################################################################################

