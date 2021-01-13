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
 Subroutine NodeMetricInterpolate(Dim,NC,NP,NConectCell,IConectCell,CellMetric,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NC,NP,I,J,K,T,Cell,ierr
 Real(8),Dimension(1:3,1:3,1:Dim)::CellMetric,NodeMetric
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:100,1:Dim)::IConectCell
 Real(8),Dimension(1:3,1:3)::U,Vt,TempMatrix2,TempMatrix1
 Real(8),Dimension(1:3)::S
 Real(8)::temp
 Integer,Dimension(1:2)::MDim1,MDim2
!*********************************************************************************************
!Part 1:
 NodeMetric(:,:,:) = 0
 Do I=1,NP
    If(NConectCell(I)==0)Cycle
    Do  J=1,NConectCell(I)
        
      !Part 2:
       Cell = IConectCell(J,I) 
       CALL SVD(CellMetric(:,:,Cell),3, S, .TRUE., U, .TRUE., Vt, ierr)
       Do K=1,3
          Do T=K,3
             temp    = U(K,T)
             U(K,T)  = U(T,K)
             U(T,K)  = temp
          EndDo
       EndDo
        
      !Part 3:
       TempMatrix2(:,:) = 0
       Do K=1,3
          If(S(K)==0.0) Cycle
          S(K) = (S(K)**(-0.5))
          TempMatrix2(K,K) = S(K)
       EndDo
        
       MDim1(:) = 3
       MDim2(:) = 3
       TempMatrix1(:,:) = 0
       Call MulMatrix(U,MDim1,TempMatrix2,MDim2,TempMatrix1)
        
       TempMatrix2(:,:) = 0
       Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
       NodeMetric(:,:,I) = NodeMetric(:,:,I)+TempMatrix2(:,:)
       
    EndDo
    
   !Part 4:
    NodeMetric(:,:,I) = (1.0/NConectCell(I))*NodeMetric(:,:,I)
    
    CALL SVD(NodeMetric(:,:,I),3, S, .TRUE., U, .TRUE., Vt, ierr)
    Do K=1,3
       Do T=K,3
          temp    = U(K,T)
          U(K,T)  = U(T,K)
          U(T,K)  = temp
       EndDo
    EndDo
    TempMatrix2(:,:) = 0
    Do K=1,3
       If(S(K)==0.0) Cycle
       S(K) = (S(K)**(-2))
       TempMatrix2(K,K) = S(K)
    EndDo
    
    MDim1(:) = 3
    MDim2(:) = 3
    TempMatrix1(:,:) = 0
    Call MulMatrix(U,MDim1,TempMatrix2,MDim2,TempMatrix1)
    
    TempMatrix2(:,:) = 0
    Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
    
    NodeMetric(:,:,I) = TempMatrix2(:,:)
 EndDo
!*********************************************************************************************
 End
!###########################################################################################