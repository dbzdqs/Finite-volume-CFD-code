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
 Subroutine ReduceMetricNodeSize(Dim,NP,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,T,K,ierr
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8),Dimension(1:2,1:2)::U,S,Vt,TempMatrix2,TempMatrix1
 Real(8),Dimension(0:2)::h
 Real(8)::temp
 Integer,Dimension(1:2)::MDim1,MDim2
 Real(8),Dimension(1:2)::EigenValue
!*********************************************************************************************
!Part 1:
 MDim1(:) = 2
 MDim2(:) = 2
 Do I=1,NP    
 
   !Part 2:
    CALL SVD(NodeMetric(:,:,I),2, EigenValue, .TRUE., U, .TRUE., Vt, ierr)
    Do K=1,2
       Do T=K,2
          temp    = U(K,T)
          U(K,T)  = U(T,K)
          U(T,K)  = temp
       EndDo
    EndDo
   
   !Part 3:
    Do J=1,2
       h(J) = 1/sqrt(EigenValue(J))
       h(J) = h(J) / 2.0
    EndDo  
   
   !Part 4:
    Do J=1,2
       EigenValue(J)    = (h(J)**-2.0)
    EndDo
    
   !Part 5:
    S(:,:) = 0
    Do J=1,2
       S(J,J)=EigenValue(J)
    EndDo
    MDim1(:) = 2
    MDim2(:) = 2
    TempMatrix1(:,:) = 0    
    Call MulMatrix(U,MDim1,S,MDim2,TempMatrix1)
    
    TempMatrix2(:,:) = 0    
    Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
    
    NodeMetric(:,:,I) = TempMatrix2(:,:)
     
 EndDo 
!*********************************************************************************************
 End
!###########################################################################################
    
    
    
    