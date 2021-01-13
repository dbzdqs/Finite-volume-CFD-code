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
!// Developed by: F. Farhadkhani, Mathmatical, Amirkabir university of Technology          //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DelCeckEBased2D(Dim,P1,P2,P3,P4,X,Y,DeL)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,P1,P2,P3,P4,X,Y
 Intent(Out  )::DeL

 Real(8)::A11,A12,A13,A21,A22,A23,A31,A32,A33,Det
 Integer::Dim,I,J,N,ME,NE,P1,P2,P3,P4,DeL
 Real(8),Dimension(1:Dim)::X,Y
!*************************************************************************************
!Part 1: 
 A11 = X(P1) - X(P4)
 A12 = Y(P1) - Y(P4)
 A13 = A11*A11 + A12*A12

 A21 = X(P2) - X(P4)
 A22 = Y(P2) - Y(P4)
 A23 = A21*A21 + A22*A22

 A31 = X(P3) - X(P4)
 A32 = Y(P3) - Y(P4)
 A33 = A31*A31 + A32*A32

 Det = A11*( A22*A33 - A23*A32 ) - A12*( A21*A33 - A23*A31  ) + A13*( A21*A32 - A22*A31 )
 
!Part 2: 
 If( Det<=0.0 )Then
 DeL=1
 Else
 DeL=-1
 Endif
!*********************************************************************************************
 End
!###########################################################################################