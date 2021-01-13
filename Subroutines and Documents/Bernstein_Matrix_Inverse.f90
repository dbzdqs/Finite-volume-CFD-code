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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Bernstein_Matrix_Inverse(N,B)
Implicit None
!********************************************************************************************* 
Intent(In   )::N
Intent(Out  )::B

Integer::N,I0,J0,N0,P,Q
Real(8)::Choose,A(1:N,1:N),B(1:N,1:N)
!*********************************************************************************************
!Part1:
Do P = 1,N
  Do Q = 1, N
   A(P,Q) = 0.0
   End Do
End Do

N0 = N - 1
!Part2;
Do J0 = 0, N0
  Do I0 = 0, J0
     A(I0+1,J0+1) = Choose ( J0, I0 ) / Choose ( N0, I0 )
  End Do
End Do
!PArt3:
Do P = 1,N
  Do Q = 1, N
   B(P,Q) = A(Q,P)
   End Do
End Do

Return
!*********************************************************************************************
End
!########################################################################################### 
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description: computes the binomial coefficient C(N,K)                                //!                                                                                      //!
!//                                                                                      //!
!// Date: August,1,2016                                                                  //!
!// Developed by: H.Moradtabrizi, Iran, Tehran, h_mtabrizi@ut.ac.ir 					 //!
!// Version: V1                                                                          //!
!// Doc ID:                                                                              //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Function choose( N, K )
 Implicit None
!********************************************************************************************* 
 Intent(In   )::N,K

 Integer::N,K,I,MN,MX
 Real(8)::Choose,Value
!*********************************************************************************************
  MN = Min ( K, N - K )

  If ( MN < 0 ) Then
    Value = 0.0D+00
  Else If ( MN == 0 ) Then
    Value = 1.0D+00
  Else
    MX = Max ( K, N - K )
    Value = Real ( MX + 1, kind = 8 )
    Do I = 2, MN
      Value = ( Value * Real ( MX + I, kind = 8 ) ) / Real ( I, kind = 8 )
    End do
  End If

  Choose = Value

 Return
!*********************************************************************************************
 End
!########################################################################################### 