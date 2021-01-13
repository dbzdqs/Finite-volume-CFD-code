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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DistributionFunction(Flag,NSegment,Lambda,Spacinge)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Flag,NSegment,Lambda
 Intent(out  )::Spacinge
 
 Integer::NSegment, Flag, Mode,I
 Real(8),Dimension(1:NSegment)::Spacinge
 Real(8)::Lambda,Sum
!*********************************************************************************************
!Part 1:
 Select case(Flag)
     
 Case(1)
  Call Distr_exp1(NSegment,Lambda, Spacinge)
 Case(2)
  Call Distr_exp2(NSegment,Lambda, Spacinge)
 Case(3)
  Call Distr_exp3(NSegment,Lambda, Spacinge)
 Case(4)
  Mode = 0
  Call Distr_cosine(MODE,NSegment,Spacinge)
 Case(5)
  Call Distr_geometric(NSegment,lambda, Spacinge)
 Case(6)
  Call Distr_uniform(NSegment, Spacinge)
 End Select
!print*,Flag,"Flag"

!!!!Part 2:
!!! Sum = 0
!!! Do I=1,NSegment+1
!!!    Sum = Sum + Spacinge(I)
!!! End Do
!!! 
!!! Do I=1,NSegment+1
!!!    Spacinge(I) = Spacinge(I)/Sum
!!! End Do
!*********************************************************************************************
 End
!###########################################################################################
