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
  Subroutine NodeMetricCorrecting(Dim,NF,NP,IDS,NConectCell,NConectEdge,IConectEdge,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NF,NP,I,J,K,T,ierr,NDirNodes,Node,Edge
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Integer,Dimension(1:Dim)::NConectCell,NConectEdge
 Real(8),Dimension(1:2,1:2)::U,Vt,CorrectedNodeDiagonalM,TempMatrix2,TempMatrix1
 Real(8),Dimension(1:2)::S,tempn1,temp1n
 Integer,Dimension(1:2)::MDim1,MDim2
 Logical::swapped
 Real(8)::temp
 Integer,Dimension(1:100)::IDirNodes
 Integer,Dimension(1:Dim,1:100)::IConectEdge
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 Do I=1,NP
    If(NConectCell(I)==0)Cycle
    CorrectedNodeDiagonalM(:,:) = 0
    
    NDirNodes = 0
    Do J=1,NConectEdge(I)
       NDirNodes = NDirNodes + 1
       Edge = IConectEdge(I,J)
       If(IDS(3,Edge)==I)Then
        Node = IDS(4,Edge)
       Else
        Node = IDS(3,Edge)
       Endif
        
       CALL SVD(NodeMetric(:,:,Node),2, S, .TRUE., U, .TRUE., Vt, ierr)
       Do K=1,2
          Do T=K,2
             temp    = U(K,T)
             U(K,T)  = U(T,K)
             U(T,K)  = temp
          EndDo
       EndDo
       
       Do k = 2-1, 1, -1
          swapped = .FALSE.
          Do t = 1, k
             If (S(t) > S(t+1)) Then
              temp    = S(t)
              S(t)    = S(t+1)
              S(t+1)  = temp
                
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
       
       TempMatrix2(:,:) = 0
       Do K=1,2
          If(S(K)==0.0) Cycle
          S(K) = (S(K)**(-0.5))
          TempMatrix2(K,K) = S(K)
       EndDo 
       CorrectedNodeDiagonalM(:,:) = CorrectedNodeDiagonalM(:,:) + TempMatrix2(:,:)
       
    EndDo
    
   !Part 3:
    CorrectedNodeDiagonalM(:,:) = (1.0/NDirNodes)*CorrectedNodeDiagonalM(:,:)
    Do K=1,2
       If(CorrectedNodeDiagonalM(K,K)==0.0) Cycle
       CorrectedNodeDiagonalM(K,K) = (CorrectedNodeDiagonalM(K,K)**(-2))
    EndDo
    
   !Part 4:
    CALL SVD(NodeMetric(:,:,I),2, S, .TRUE., U, .TRUE., Vt, ierr)
    Do K=1,2
       Do T=K,2
          temp    = U(K,T)
          U(K,T)  = U(T,K)
          U(T,K)  = temp
       EndDo
    EndDo
    
    Do k = 2-1, 1, -1
       swapped = .FALSE.
       Do t = 1, k
       
          If (S(t) > S(t+1)) Then
           temp    = S(t)
           S(t)    = S(t+1)
           S(t+1)  = temp
                
           tempn1(:)    = U(t,:)
           U(t,:)    = U(t+1,:)
           U(t+1,:)  = tempn1(:)
                
           temp1n(:)    = Vt(:,t)
           Vt(:,t)    = Vt(:,t+1)
           Vt(:,t+1)  = temp1n(:)
                
           swapped = .TRUE.
          EndIf
       
       EndDo
       IF (.NOT. swapped) Exit
    EndDo
    
   !Part 5:
    MDim1(:) = 2
    MDim2(:) = 2
    TempMatrix1(:,:) = 0
    Call MulMatrix(U,MDim1,CorrectedNodeDiagonalM,MDim2,TempMatrix1)
    
    TempMatrix2(:,:) = 0
    Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
    NodeMetric(:,:,I) = TempMatrix2(:,:)
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################