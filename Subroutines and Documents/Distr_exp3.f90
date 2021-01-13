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
!// Developed by: M. Gharibi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
SUBROUTINE Distr_exp3(N,Lambda,Spacinge)
 Implicit None
!*********************************************************************************************
 INTEGER::N,M,I,J
 Real*8,Dimension(1:N+1)::Spacinge,X
 Real*8::Lambda
 INTRINSIC::EXP
!*********************************************************************************************
!part1
 M=Float(N/2)
 Do I=1,M
    X(I)=I-1
 End Do
 
!part2
 Do J=1,M
    Spacinge(J)= 0.5*(EXP(Lambda*(X(J)-M)))
 End Do
 
!part3
 Do I=M,N
    X(I)=I
 End Do
 
!part4
 Do J=M,N
    Spacinge(J+1)= 1-0.5*(EXP(-Lambda*(X(J)-M)))
 End Do
!*********************************************************************************************
 End
!###########################################################################################