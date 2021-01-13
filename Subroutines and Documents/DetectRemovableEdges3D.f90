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
 Subroutine DetectRemovableEdges3D(Dim,CF,NP,NF,NC,IDS,FaceType,NFR,BeginOfReg,X,Y,Z,Bound,NConectCell,IConectCell,NFace_Cell,IFace_Cell,NContractedEdges,IContractedEdges)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,NApexVertices,NOrderedVertices,I,T,K,Vertice,NContractedEdges,temp
 Integer,Dimension(1:100)::NFR,BeginOfReg
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,NFace_Cell,LayerIndex,NConnectedPoints,IApexVertices,NConectCell,IOrderedVertices,VST !Vertice Status: 0:unknown - 1:included - 2:excluded
 Real(8),Dimension(1:3,1:3,1:Dim)::CellMetric,NodeMetric,CoarsedNodeMetric
 Integer,PARAMETER:: Unknown = 0, Included=1, Excluded=2
 Real(8),Dimension(1:Dim)::X,Y,Z,VNN
 Logical,Dimension(1:Dim)::Bound
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints,IConectCell
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Logical::swapped
 Real(8)::CF
!*********************************************************************************************
 
!Part 1:
 Call FindConnectedEdges3D(Dim,NF,NP,IDS,FaceType,NConnectedPoints,IConnectedPoints)
 
!Part 2:
 Call MetricDefine3D(Dim,NC,NF,NP,IDS,FaceType,X,Y,Z,NFace_Cell,IFace_Cell,NConectCell,IConectCell,NConnectedPoints,IConnectedPoints,CellMetric,NodeMetric)
 
!Part 3:
 Call DefineBoundaryLayerIndexes3D(Dim,NP,NF,BeginOfReg,NFR,NC,IDS,FaceType,X,Y,Z,Bound,NConnectedPoints,IConnectedPoints,NodeMetric,LayerIndex)
 
!Part 4: 
 NOrderedVertices = 0
 Do I=1,NP
    If(NConnectedPoints(I)==0 .OR. LayerIndex(I)<=0)Cycle
    NOrderedVertices = NOrderedVertices + 1
    IOrderedVertices(NOrderedVertices)  = I    
 EndDo
 
!Part 5:
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
 
!Part 6:
 Do I=1,NP
    If(NConnectedPoints(I)==0 .OR. LayerIndex(I)>0)Cycle
    NOrderedVertices = NOrderedVertices + 1
    IOrderedVertices(NOrderedVertices)  = I    
 EndDo
 
!Part 7:
 Call CalcVerticesNN3D(Dim,NP,X,Y,Z,NConnectedPoints,IConnectedPoints,VNN)
 
!Part 8:
 Call NodeMetricCoarsening3D(Dim,NP,CF,NConectCell,NodeMetric,CoarsedNodeMetric)
 
!Part 9:
 Call HShockCorrection3D(Dim,NF,NP,IDS,FaceType,X,Y,Z,CoarsedNodeMetric)
 
!Part 10:
 Call ReduceMetricNodeSize3D(Dim,NP,CoarsedNodeMetric)
 
!Part 11:
 VST(:) = Unknown
 Call FindBoundApexVertices(Dim,NF,NP,IDS,FaceType,NApexVertices,IApexVertices,Bound,X,Y,Z)
 Do I=1,NApexVertices
     VST(IApexVertices(I))=Included
 EndDo
 
!Part 12:
 NContractedEdges = 0
 Do I=1,NOrderedVertices
    Vertice = IOrderedVertices(I) 
    
   !Part 13:
    If(VST(Vertice)==Excluded ) Cycle
    If(LayerIndex(Vertice)>0) VST(Vertice)=Included
    Call CheckNeighborVertices3D(Dim,NP,NF,NC,IDS,FaceType,NConectCell,IConectCell,NFace_Cell,IFace_Cell,Bound,Vertice,NConnectedPoints,IConnectedPoints,VNN,NodeMetric,CoarsedNodeMetric,X,Y,Z,LayerIndex,VST,NContractedEdges,IContractedEdges)
    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################