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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine KFi_Dif3d(Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DFiX_F,DFiY_F,DFiZ_F,Sigk,SigFi,&
                        MR,Mu,Mut,Dift)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NFW1,NFW2,NF,NF1,NF2,IDS,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DFiX_F,DFiY_F,DFiZ_F,Sigk,SigFi,MR,Mu,Mut
 Intent(Out  )::Dift

 Integer::Dim,I,NC,NFW1,NFW2,NF,NF1,NF2,ME,NE
 Real(8)::Mu_k,Mu_w,Sigk,SigFi,MR,F1,F2
 Real(8),Dimension(1:Dim)::Mu,Mut,NX,NY,NZ,DKX_F,DKY_F,DKZ_F,DFiX_F,DFiY_F,DFiZ_F
 Real(8),Dimension(1:2,1:Dim)::Dift
 Integer,Dimension(1:6,1:Dim)::IDS
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
    Mu_k = 0.5*( Mu(ME) + Mut(ME)*Sigk  + Mu(NE) + Mut(NE)*Sigk  )
	Mu_w = 0.5*( Mu(ME) + Mut(ME)*SigFi + Mu(NE) + Mut(NE)*SigFi )
    
	F1 = Mu_k * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) + DKZ_F(I)  *NZ(I) ) 
	F2 = Mu_w * ( DFiX_F(I) *NX(I) + DFiY_F(I) *NY(I) + DFiZ_F(I) *NZ(I) ) 

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
	F1 = Mu(ME) * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) + DKZ_F(I)  *NZ(I) ) 
	F2 = Mu(ME) * ( DFiX_F(I) *NX(I) + DFiY_F(I) *NY(I) + DFiZ_F(I) *NZ(I) ) 

   !Part 9:
    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
    
 End do
 
!Part 10:
 Do I=NFW2+1,NF

    ME = IDS(1,I)
  
    Mu_k = Mu(ME) + Mut(ME)*Sigk
	Mu_w = Mu(ME) + Mut(ME)*SigFi
    
	F1 = Mu_k * ( DKX_F(I)  *NX(I) + DKY_F(I)  *NY(I) + DKZ_F(I)  *NZ(I) ) 
	F2 = Mu_w * ( DFiX_F(I) *NX(I) + DFiY_F(I) *NY(I) + DFiZ_F(I) *NZ(I) ) 

    Dift(1,ME) = Dift(1,ME) + F1
    Dift(2,ME) = Dift(2,ME) + F2
 
 End do

!Part 11:
 Do I=1,NC
    Dift(1,I) = -MR*Dift(1,I)
	Dift(2,I) = -MR*Dift(2,I)
 End do
!*********************************************************************************************
 End
!###########################################################################################

