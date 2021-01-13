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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine LES_DSmag_IsoSGS3D(Dim,NC,NF1,NF2,IDS,MR,Vol,Rho,Rhohat,Sabs,Sabshat,Lkk,Taukk)
 
 Implicit None
!********************************************************************************************* 
 Intent(IN)::Dim,NC,NF1,NF2,MR,Vol,IDS,Rho,Rhohat,Sabs,Sabshat,Lkk
 Intent(INOUT)::Taukk

 Integer::Dim,NC,NF1,NF2,I,Allocatestatus,DeAllocatestatus,ME
 Real(8)::MR,Alpha,CI

 Real(8),Dimension(1:Dim)::Vol,Taukk,Rho,Rhohat,Sabs,Sabshat,Lkk
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::Beta,Betahat,A_Bhat,Lkkhat,ABhat
!********************************************************************************************* 
!Part 1:
 Do I=1,NC
    Beta(I) = 2*Rho(I) * (Vol(I)**0.6666) * Sabs(I)*Sabs(I)
 End Do

!Part 2:
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Beta,Vol,Betahat)
      
!Part 3:
 Do I=1,NC
    Alpha = 2*Rhohat(I)* 4*(Vol(I)**0.6666) * Sabshat(I)*Sabshat(I)

    A_Bhat(I) = Alpha-Betahat(I)
 End Do

!Part 4:
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,A_Bhat,Vol,ABhat)
 call LES_Explicit_Averaging3D(Dim,NC,NF1,NF2,IDS,Lkk,Vol,Lkkhat)

 Do I=1,NC

   !Part 5:
    CI = Lkkhat(I)/ABhat(I)

    If (CI<0.0)Then
      CI=0.0;
    Else If (CI>0.0066) Then
      CI=0.0066
    Else If (CI>=0.0 .And. CI<=0.0066) Then
      CI=CI
    Else 
      CI=0.0066
    End If

   !Part 6:
    Taukk(I) = (1.0/MR) * CI * 2*Rho(I)* (Vol(I)**0.6666) * Sabs(I)*Sabs(I)

 End Do
!*********************************************************************************************
 End 
!########################################################################################### 
