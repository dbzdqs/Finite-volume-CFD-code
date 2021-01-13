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
!// Date: May., 15, 2016                                                                   //!
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
  Subroutine DetectRemovableEdges(Dim,CF,NP,NF,NC,IDS,X,Y,NEdgeOfCell,IEdgeOfCell,NContractedEdges,IContractedEdges)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,NOrderedVertices,I,K,T,Vertice,NContractedEdges,temp
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,VNN,ElipsoidRatio
 Logical,Dimension(1:Dim)::Bound,StretchNode
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:Dim)::NConnectedEdges,NConnectedNodes,IOrderedVertices,LayerIndex,VST,NConectCell
 Integer,Dimension(1:50,1:Dim)::IConectCell
 Integer,Dimension(1:4,1:Dim)::IEdgeOfCell,IPoint_Cell
 Integer,Dimension(1:Dim)::NEdgeOfCell,NPoint_Cell
 Integer,PARAMETER:: Unknown = 0, Included=1, Excluded=2
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
 Logical::swapped
 Real(8)::CF
 Real(8),Dimension(1:2,1:2,1:Dim)::CellMetric
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric,CoarsedNodeMetric
!*********************************************************************************************
!Part 1:
 Call PointOfCellV2(Dim,NF,NC,IDS,NEdgeOfCell,IEdgeOfCell,NPoint_Cell,IPoint_Cell)
!Part 2:
 Call ConectedCellV2(Dim,NC,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NP)
!Part 3:
 Call ConectedEdgeOfPoint(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges)
!Part 4:
 Call FindConnectedNodes(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges,NConnectedNodes,IConnectedNodes) 
!Part 5:
 Call BoundPointLabelingV2(Dim,Bound,NF,IDS)
!Part 6:
 Call MetricDefine(Dim,NC,NF,NP,IDS,X,Y,NEdgeOfCell,IEdgeOfCell,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,CellMetric,NodeMetric)
!Part 7:
 Call NodeMetricCorrecting(Dim,NF,NP,IDS,NConectCell,NConnectedEdges,IConnectedEdges,NodeMetric)
!Part 8:
 Call DefineBoundaryLayerIndexes(Dim,NP,NF,NC,IDS,X,Y,Bound,NConnectedEdges,IConnectedEdges,NConnectedNodes,IConnectedNodes,NodeMetric,LayerIndex)
 
!Part 9: 
 NOrderedVertices = 0
 Do I=1,NP
    If(NConnectedEdges(I)==0 .OR. LayerIndex(I)<=0)Cycle
    NOrderedVertices = NOrderedVertices + 1
    IOrderedVertices(NOrderedVertices)  = I    
 EndDo
 
!Part 10:
 Do k = NOrderedVertices-1, 1, -1
    swapped = .FALSE.
    Do t = 1, k
       
       If (abs(LayerIndex(IOrderedVertices(t))) > abs(LayerIndex(IOrderedVertices(t+1)))) Then
        temp    = IOrderedVertices(t)
        IOrderedVertices(t)    = IOrderedVertices(t+1)
        IOrderedVertices(t+1)  = temp
        swapped = .TRUE.
       EndIf
       
    EndDo
    IF (.NOT. swapped) Exit
 EndDo
 
!Part 11:
 Do I=1,NP
    If(NConnectedEdges(I)==0 .OR. LayerIndex(I)>0)Cycle
    NOrderedVertices = NOrderedVertices + 1
    IOrderedVertices(NOrderedVertices)  = I    
 EndDo

!Part 12:
 Call CalcVerticesNN(Dim,NP,X,Y,NConnectedNodes,IConnectedNodes,VNN)
 
!Part 13:
 Call NodeMetricCoarsening(Dim,NP,CF,NConectCell,NodeMetric,CoarsedNodeMetric)
 
!Part 14:
 Call HShockCorrection(Dim,NF,NP,IDS,X,Y,NodeMetric) 
 
!Part 15:
 Call ReduceMetricNodeSize(Dim,NP,CoarsedNodeMetric)
 
!Part 16:
 NContractedEdges = 0
 VST(:) = Unknown
 Do I=1,NOrderedVertices
    
   !Part 17:
    Vertice = IOrderedVertices(I) 
    If(VST(Vertice)==Excluded ) Cycle
    Call CheckNeighborVertices(Dim,NP,NF,NC,Bound,Vertice,NConnectedNodes,IConnectedNodes,VNN,NodeMetric(:,:,1:NP),CoarsedNodeMetric(:,:,1:NP),X,Y,LayerIndex,VST,NContractedEdges,IContractedEdges)
    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################