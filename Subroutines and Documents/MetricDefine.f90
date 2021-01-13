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
  Subroutine MetricDefine(Dim,NC,NF,NP,IDS,X,Y,NFace_Cell,IFace_Cell,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,CellMetric,NodeMetric)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NC,NF,NP,NTRemNode,RNode,NTempCellEdges,NCellEdges,J2,TempNode,T,J,K,P1,P2,I,Cell,ierr,NCell
 Integer,Dimension(1:2,1:4)::ITempCellEdges,ICellEdges
 Real(8),Dimension(1:Dim)::X,Y
 Integer,Dimension(1:4,1:Dim)::IFace_Cell,IPoint_Cell
 Integer,Dimension(1:Dim)::NFace_Cell,NPoint_Cell
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:50,1:Dim)::IConectCell
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:2,1:2,1:Dim)::CellMetric
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8)::temp
 Logical::exists
 Integer,Dimension(1:2)::MDim1,MDim2,TRemNode
 Real(8),Dimension(1:2,1)::TNodeCP1P2
 Real(8),Dimension(1,1:2)::NodeCP1P2
 Real(8),Dimension(1:2)::S
 Real(8),Dimension(1:2,1:2)::U,Vt,TempMatrix1,TempMatrix,TempMatrix2,TempCellMetric
!*********************************************************************************************
!Part 1:
 NodeMetric(:,:,:) = 0
 Do I=1,NP
    If(NConectCell(I)==0)Cycle

   !Part 2:
    NCell = 0
    Do  J=1,NConectCell(I)
       Cell = IConectCell(J,I) 
      
      !Part 3:
       NTRemNode   = 0
       TRemNode(:) = 0
       If(NFace_Cell(Cell)==4)Then
           
        Do K=1,NFace_Cell(Cell)
           If(IDS(3,IFace_Cell(K,Cell))==I)Then
            NTRemNode = NTRemNode + 1
            TRemNode(NTRemNode) = IDS(4,IFace_Cell(K,Cell))
           ElseIf(IDS(4,IFace_Cell(K,Cell))==I)Then
            NTRemNode = NTRemNode + 1
            TRemNode(NTRemNode) = IDS(3,IFace_Cell(K,Cell))
           Endif
        EndDo
       Else
        NTRemNode = 1
       Endif
       
      !Part 4:
       Do K=1,NTRemNode
          RNode = TRemNode(K)
          NTempCellEdges = NFace_Cell(Cell)
           
          TempNode = 0
          Do T=1,NTempCellEdges
             If(IDS(3,IFace_Cell(T,Cell))==RNode .AND. IDS(4,IFace_Cell(T,Cell))/=I )Then
              TempNode = IDS(4,IFace_Cell(T,Cell))
             ElseIf(IDS(4,IFace_Cell(T,Cell))==RNode .AND. IDS(3,IFace_Cell(T,Cell))/=I )Then
              TempNode = IDS(3,IFace_Cell(T,Cell))
             EndIf
          EndDo
        
         !Part 5:
          Do T=1,NFace_Cell(Cell)
             If(IDS(1,IFace_Cell(T,Cell))==Cell)Then
              P1 = 3
              P2 = 4
             Else
              P1 = 4
              P2 = 3
             Endif
               
             If(IDS(P1,IFace_Cell(T,Cell))==RNode)Then
              ITempCellEdges(1,T) = TempNode
             Else
              ITempCellEdges(1,T) = IDS(P1,IFace_Cell(T,Cell))
             Endif
               
             If(IDS(P2,IFace_Cell(T,Cell))==RNode)Then
              ITempCellEdges(2,T) = TempNode
             Else
              ITempCellEdges(2,T) = IDS(P2,IFace_Cell(T,Cell))
             Endif
          EndDo
          
         !Part 6:
          NCellEdges = 0
          Do T=1,NTempCellEdges
             exists = .FALSE.
             Do J2=1,NCellEdges
                If( ITempCellEdges(2,T)==ITempCellEdges(1,T) )Then
                 exists = .TRUE.
                 exit
                Endif
             EndDo
             If(.NOT. exists)Then
              NCellEdges = NCellEdges + 1
              ICellEdges(1,NCellEdges) = ITempCellEdges(1,T)
              ICellEdges(2,NCellEdges) = ITempCellEdges(2,T)
             Endif
          EndDo
           
         !Part 7:
          TempCellMetric(:,:) = 0
          Do T=1,NCellEdges
             P1 = ICellEdges(1,T)
             P2 = ICellEdges(2,T)
               
             NodeCP1P2(:,:)  = 0
             TNodeCP1P2(:,:) = 0
             NodeCP1P2(1,1)  = X(P2)-X(P1)
             NodeCP1P2(1,2)  = Y(P2)-Y(P1)
             TNodeCP1P2(1,1) = NodeCP1P2(1,1)
             TNodeCP1P2(2,1) = NodeCP1P2(1,2)
        
             MDim1(1) = 2
             MDim1(2) = 1     
             MDim2(1) = 1
             MDim2(2) = 2     
             TempMatrix(:,:) = 0
             Call MulMatrix(NodeCP1P2,MDim1,TNodeCP1P2,MDim2,TempMatrix)
       
             TempCellMetric(:,:) = TempCellMetric(:,:) + TempMatrix(:,:)  
          EndDo
           
         !Part 8:
          TempMatrix(:,:) = 0
          Call CalcMatrixInverseV2(TempCellMetric(:,:),TempMatrix(:,:),2)
          TempCellMetric(:,:)     = TempMatrix(:,:)
          TempCellMetric(1:2,1:2) = ((1.5)*TempCellMetric(1:2,1:2))
           
           
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!Node Metric Interpolate !!!! 
         !Part 9:
          CALL SVD(TempCellMetric(:,:),2, S, .TRUE., U, .TRUE., Vt, ierr)
          Do J2=1,2
             Do T=J2,2
                temp     = U(J2,T)
                U(J2,T)  = U(T,J2)
                U(T,J2)  = temp
             EndDo
          EndDo
        
         !Part 10:
          TempMatrix2(:,:) = 0
          Do J2=1,2
             If(S(J2)==0.0) Cycle
             S(J2) = (S(J2)**(-0.5))
             TempMatrix2(J2,J2) = S(J2)
          EndDo
        
          MDim1(:) = 2
          MDim2(:) = 2
          TempMatrix1(:,:) = 0
          Call MulMatrix(U,MDim1,TempMatrix2,MDim2,TempMatrix1)
        
          TempMatrix2(:,:) = 0
          Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
       
          NCell = NCell + 1
       
          NodeMetric(:,:,I) = NodeMetric(:,:,I) + TempMatrix2(:,:)
          !!!!Node Metric Interpolate !!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
       EndDo
       
    EndDo
   
   !Part 11:
    NodeMetric(:,:,I) = (1.0/NCell)*NodeMetric(:,:,I)
   
    CALL SVD(NodeMetric(:,:,I),2, S, .TRUE., U, .TRUE., Vt, ierr)
    Do K=1,2
       Do T=K,2
          temp    = U(K,T)
          U(K,T)  = U(T,K)
          U(T,K)  = temp
       EndDo
    EndDo
    
    TempMatrix2(:,:) = 0
    Do K=1,2
       If(S(K)==0.0) Cycle
       S(K) = (S(K)**(-2))
       TempMatrix2(K,K) = S(K)
    EndDo
    
    MDim1(:) = 2
    MDim2(:) = 2
    TempMatrix1(:,:) = 0
    Call MulMatrix(U,MDim1,TempMatrix2,MDim2,TempMatrix1)
    
    TempMatrix2(:,:) = 0
    Call MulMatrix(TempMatrix1,MDim1,Vt,MDim2,TempMatrix2)
    
    NodeMetric(:,:,I) = TempMatrix2(:,:)
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################