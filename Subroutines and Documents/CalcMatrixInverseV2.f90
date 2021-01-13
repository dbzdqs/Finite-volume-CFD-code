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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CalcMatrixInverseV2(a,c,n)
!*********************************************************************************************
 implicit none 
 integer n
 double precision a(n,n), c(n,n)
 Real(16):: L(n,n), U(n,n), b(n), d(n), x(n)
 Real(16)::coeff
 integer i, j, k
!*********************************************************************************************
 L = 0.0
 U = 0.0
 b = 0.0
 
 Do K=1,n-1
    Do I=K+1,n
       coeff  = a(I,K)/a(K,K)
       L(I,K) = coeff
       Do j=k+1,n
          a(I,J) = a(I,J)-coeff*a(K,J)
       EndDo
    EndDo
 EndDo
 
 Do I=1,n
    L(I,I) = 1.0
 EndDo
 
 Do J=1,n
    Do I=1,J
       U(I,J) = a(I,J)
    EndDo
 EndDo
 
 Do K=1,n
    b(K) = 1.0
    d(1) = b(1)
    
    Do I=2,n
       d(I) = b(I)
       Do J=1,I-1
          d(I) = d(I) - L(I,J)*d(J)
       EndDo
    EndDo
    
    x(n) = d(n)/U(n,n)
    Do I = n-1,1,-1
       x(I) = d(I)
       Do J=n,I+1,-1
          x(I) = x(I)-U(I,J)*x(J)
       EndDo
       x(I) = x(I)/u(I,I)
    EndDo
    
    Do I=1,n
       c(I,K) = x(I)
    EndDo
    
    b(K)=0.0
 EndDo
!*********************************************************************************************
 End
!###########################################################################################