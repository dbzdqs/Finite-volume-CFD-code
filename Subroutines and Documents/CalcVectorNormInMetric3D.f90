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
 Subroutine CalcVectorNormInMetric3D(Metric,Vector,Norm)
 Implicit None
!*********************************************************************************************
 Real(8),Dimension(1:3,1:3)::Metric
 Real(8)::Norm
 Real(8),Dimension(1,3)::Vector
 Real(8),Dimension(3,1)::VectorT
 Integer,Dimension(1:3)::MDim1,MDim2
 Real(8),Dimension(1)::TempMatrix4
 Real(8),Dimension(1:3,1)::TempMatrix2
!*********************************************************************************************
!Part 1:
 VectorT(1,1) = Vector(1,1)
 VectorT(2,1) = Vector(1,2)
 VectorT(3,1) = Vector(1,3)
        
!Part 2:
 MDim1(1) = 1
 MDim1(2) = 3
 MDim2(1) = 3
 MDim2(2) = 3
 TempMatrix2(1:3,1) = 0
 Call MulMatrix(VectorT,MDim1,Metric(:,:),MDim2,TempMatrix2)

!Part 3:
 MDim1(1) = 1
 MDim1(2) = 3
 MDim2(1) = 3
 MDim2(2) = 1
 TempMatrix4(1) = 0
 Call MulMatrix(TempMatrix2,MDim1,Vector,MDim2,TempMatrix4)
        
!Part 4:
 Norm = sqrt(TempMatrix4(1))
 
!*********************************************************************************************
End    
