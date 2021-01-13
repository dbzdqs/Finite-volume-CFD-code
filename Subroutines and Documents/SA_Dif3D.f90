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
 Subroutine SA_Dif3D(Dim,NC,NF1,NF2,NF,MR,Cb2,Sigma,NX,NY,NZ,Vol,WTNP1,WNP1,Mu,WTB,IDS,&
                     DNuX,DNuY,DNuZ,DRNuX,DRNuY,DRNuZ,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc,Dift)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,MR,Cb2,Sigma,NX,NY,NZ,Vol,WTNP1,WNP1,Mu,WTB,IDS,&
                DNuX,DNuY,DNuZ,DRNuX,DRNuY,DRNuZ,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc
 Intent(Out  )::Dift

 Integer::Dim,I,NC,NF1,NF2,NF,ME,NE
 Real(8)::T1,Mu_eff,Del_Nu,MR,Cb2,Sigma
 Real(8),Dimension(1:Dim)::NX,NY,NZ,Vol,Mu,DNuX,DNuY,DNuZ,DRNuX,DRNuY,DRNuZ,DNuXc,DNuYc,DNuZc,DRNuXc,DRNuYc,DRNuZc
 Real(8),Dimension(1:2,1:Dim)::WTNP1,WTB,Dift
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::WNP1
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    Dift(1,I) = Vol(I) * Cb2 * ( DNuXc(I)*DRNuXc(I) + DNuYc(I)*DRNuYc(I) + DNuZc(I)*DRNuZc(I) )
 End do

!Part 2:
 Do I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    Mu_eff = Mu(ME) + WTB(1,I)
    T1 = Mu_eff * ( DNuX(I)*NX(I) + DNuY(I)*NY(I) + DNuZ(I)*NZ(I) )

    Dift(1,ME) = Dift(1,ME) + T1

 End do

!Part 5:
 Do I=NF1+1,NF2

   !Part 6:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 7:
    Mu_eff = 0.5*( Mu(ME)+WTNP1(1,ME) + Mu(NE)+WTNP1(1,NE) )

    T1 = Mu_eff * ( DNuX(I)*NX(I) + DNuY(I)*NY(I) + DNuZ(I)*NZ(I) )

   !Part 8:
    Dift(1,ME) = Dift(1,ME) + T1
    Dift(1,NE) = Dift(1,NE) - T1

 End do

!Part 9:
 Do I=1,NC
    Dift(1,I) = MR * Dift(1,I) / Sigma
 End do
!**********************************************************************************************
 End
!##############################################################################################


