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
 SUBROUTINE Distr_cosine(MODE,N,Spacinge)
 Implicit None
!*********************************************************************************************
 INTEGER::N,N1,I,J,K,B,MODE
 Real*8,Dimension(1:N+1)::Spacinge,X
 Real*8::Pi,A
 INTRINSIC::DSIN
!*********************************************************************************************
!part1
 N1=N-1
 A=1
 B=0
 Pi= 3.1416
 
 !!!!!!!!!!arbitary cosine distribution
!part2
 If(MODE.EQ.0)Then
  Do I=1,N+1
     X(I)=(B-(0.5*A))+((A/N1)*(I-1))
  End Do
  
  Do J=1,N
     Spacinge(J)= 0.5*(1+DSIN(((X(J)-B)/A)*Pi))
  End Do
  
  Spacinge(1)=0
  
!!!!!!!!!!Raised cosine distribution
!part3
 Else If(MODE.EQ.1)Then
  Do I=1,N+1
     X(I)=(B-A)+(2*A/N1)*(I-1)
  End Do
  
  Do J=1,N
     Spacinge(J+1)= 0.5*(1+((X(J)-B)/A)+(1/Pi)*DSIN(((X(J)-B)/A)*Pi))
  End Do
  Spacinge(1)=0
  
 End If
!*********************************************************************************************
 End
!###########################################################################################