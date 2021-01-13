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
 Subroutine CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric,Node1,Node2,h)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,Node1,Node2
 Real(8),Dimension(1:2,1:2)::NodeMetric
 Real(8)::h,Norm
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1,2)::UnitV
 Real(8),Dimension(2,1)::UnitVt
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1)::TempMatrix1,TempMatrix4
 Real(8),Dimension(1:2,1)::TempMatrix2
 Real(8),Dimension(1:2,1:10)::PVector
!*********************************************************************************************
!Part 1:
 PVector(1,1)  = X(Node2)-X(Node1)
 PVector(2,1)  = Y(Node2)-Y(Node1)
 Call CalcVectorNormInMetric(NodeMetric(:,:),PVector(:,1),Norm)
 
 If(Norm==0)Then
  h = 0.0
  return
 Endif
        
!Part 2:
 UnitV(1,1)  = (X(Node2)-X(Node1))/Norm
 UnitV(1,2)  = (Y(Node2)-Y(Node1))/Norm
 UnitVt(1,1) = UnitV(1,1)
 UnitVt(2,1) = UnitV(1,2)
        
!Part 3:
 MDim1(1) = 1
 MDim1(2) = 2
 MDim2(1) = 2
 MDim2(2) = 1
 TempMatrix1(1) = 0
 Call MulMatrix(UnitVt,MDim1,UnitV,MDim2,TempMatrix1)
        
 MDim1(1) = 1
 MDim1(2) = 2
 MDim2(1) = 2
 MDim2(2) = 2
 TempMatrix2(1:2,1) = 0
 Call MulMatrix(UnitVt,MDim1,NodeMetric(:,:),MDim2,TempMatrix2)
            
 MDim1(1) = 1
 MDim1(2) = 2
 MDim2(1) = 2
 MDim2(2) = 1
 TempMatrix4(1) = 0
 Call MulMatrix(TempMatrix2,MDim1,UnitV,MDim2,TempMatrix4)
            
 h = sqrt(TempMatrix1(1))/sqrt(TempMatrix4(1))
 
!*********************************************************************************************
End
!###########################################################################################