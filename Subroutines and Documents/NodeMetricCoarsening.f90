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
 Subroutine NodeMetricCoarsening(Dim,NP,CF,NConectCell,NodeMetric,CoarsedNodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,J3,T,K,ierr
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric,CoarsedNodeMetric
 Real(8),Dimension(1:2,1:2)::U,S,Vt,TempMatrix2,TempMatrix1
 Real(8),Dimension(0:2)::h
 Real(8)::CF,temp
 Integer,Dimension(1:2)::MDim1,MDim2
 Integer,Dimension(1:Dim)::NConectCell
 Real(8),Dimension(1:2)::EigenValue,tempn1,temp1n
 Logical::swapped
!*********************************************************************************************
!Part 1:
 MDim1(:) = 2
 MDim2(:) = 2
 Do I=1,NP    
    If(NConectCell(I)==0)Cycle
 
   !Part 2:
    CALL SVD(NodeMetric(:,:,I),2, EigenValue, .TRUE., U, .TRUE., Vt, ierr)
    Do K=1,2
       Do T=K,2
          temp    = U(K,T)
          U(K,T)  = U(T,K)
          U(T,K)  = temp
       EndDo
    EndDo
    Do j3 = 2-1, 1, -1
       swapped = .FALSE.
       Do t = 1, j3
          If (EigenValue(t) < EigenValue(t+1)) Then
           temp             = EigenValue(t)
           EigenValue(t)    = EigenValue(t+1)
           EigenValue(t+1)  = temp
         
           tempn1(:)    = U(t,:)
           U(t,:)       = U(t+1,:)
           U(t+1,:)     = tempn1(:)
         
           temp1n(:)    = Vt(:,t)
           Vt(:,t)      = Vt(:,t+1)
           Vt(:,t+1)    = temp1n(:)
         
           swapped = .TRUE.
          EndIf
       EndDo
       IF (.NOT. swapped) Exit
    EndDo
     
    Do J=1,2
       h(J) = (EigenValue(J)**(-0.5))
    EndDo
    h(0) = CF*h(1)
     
   !Part 3:
    Do J=1,2
       h(J) = MAX(h(J),MIN( (CF*h(J)) , h(J-1) ))
       EigenValue(J)    = (h(J)**-2.0)
    EndDo
    
   !Part 4:
    S(:,:) = 0
    Do J=1,2
       S(J,J)=EigenValue(J)
    EndDo
    
   !Part 5:
    MDim1(:) = 2
    MDim2(:) = 2
    TempMatrix1(:,:) = 0    
    Call MulMatrix(U,MDim1,S,MDim2,TempMatrix1)
    
    TempMatrix2(:,:) = 0    
    Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
    
    CoarsedNodeMetric(:,:,I) = TempMatrix2(:,:)
     
 EndDo 
!*********************************************************************************************
 End
!###########################################################################################