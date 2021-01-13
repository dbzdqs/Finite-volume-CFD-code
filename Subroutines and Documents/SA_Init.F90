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
 Subroutine SA_Init(Dim,NC,Mu0,R0,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mut0,Nuhat0,WTNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,Mu0,R0
 Intent(Out  )::Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mut0,Nuhat0,WTNP1,Mut

 Integer::Dim,I,J,NC
 Real(8)::Chi,Nut0,Cb1,Cb2,Kei,Sigma,Cw1,Cw2,Cw3,Cv1,Mu0,R0,Xj,Yj,Xi,Yi,Dmin,Dis,DX,DY,Nu0,&
          Nuhat0,Mut0,Cv13,Chi3,Fv1
 Real(8),Dimension(1:Dim)::WTNP1,Mut
!*********************************************************************************************	
!Part 1:
 Cb1   = 0.1355
 Cb2   = 0.622
 Kei   = 0.41
 Sigma = 2/3.
 Cw1   = cb1/(Kei*Kei) + (1+cb2)/sigma
 Cw2   = 0.3
 Cw3   = 2.
 Cv1   = 7.1

!Part 2:
 Chi = 3.

!Part 3:
 Nuhat0 = Chi*R0*Mu0

!Part 4:
 Cv13 = Cv1*Cv1*Cv1
 Chi3 = Chi*Chi*Chi
 Fv1  = Chi3/(Chi3+Cv13)  

!Part 5:
 Mut0 = Fv1*R0*Nuhat0 

!Part 6:
 Do I=1,NC
    WTNP1(I) = Nuhat0
	Mut(I)   = Mut0
 End do

!*********************************************************************************************
 End
!###########################################################################################
