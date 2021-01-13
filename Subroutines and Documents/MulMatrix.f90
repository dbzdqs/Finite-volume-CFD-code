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
!// Date: Feb., 10, 2015                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine MulMatrix(Matrix1,MDim1,Matrix2,MDim2,ResultMatrix)
 Implicit None
!*********************************************************************************************
 Intent(IN)::Matrix1,MDim1,Matrix2,MDim2
 Intent(INOUT)::ResultMatrix

 Integer::I,J,K
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1:MDim1(2),1:MDim1(1))::Matrix1
 Real(8),Dimension(1:MDim2(2),1:MDim2(1))::Matrix2
 Real(8),Dimension(1:MDim2(2),1:MDim1(1))::ResultMatrix
!*********************************************************************************************

!Part 1:
 If(MDim1(2)/=MDim2(1))Then
  Print *,'Invalid Dimensions'
  pause
  stop
 Endif
  
!Part 2:
 ResultMatrix(:,:) = 0
 
!Part 3:
 Do I=1,MDim1(1)
    Do J=1,MDim2(2)
       Do K=1,MDim1(2)
          ResultMatrix(J,I) = ResultMatrix(J,I)+(Matrix1(K,I)*Matrix2(J,K))
       EndDo
    EndDo
 EndDo
 
!*********************************************************************************************
 End   
!###########################################################################################