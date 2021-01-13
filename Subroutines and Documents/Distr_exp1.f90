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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: M. Gharibi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 SUBROUTINE Distr_exp1(NSegment,Lambda,Spacinge)
 Implicit None
!*********************************************************************************************
 INTEGER::NSegment,J
 Real(8),Dimension(1:NSegment)::Spacinge
 Real(8)::Lambda
!*********************************************************************************************
!!part1
! Do I=1,N+1
!    X(I)=I-1
! End Do
!
!!part2
! Do J=1,N+1
!    Spacinge(J)= 1-(EXP(-Lambda*X(J)))
! End Do
 Lambda=0.1
 Do J=1,NSegment
    Spacinge(J)= 1.0-exp( -Lambda*(J-1) )
    !print*,Spacinge(J),j
 End Do
 !pause
!*********************************************************************************************
 End
!###########################################################################################