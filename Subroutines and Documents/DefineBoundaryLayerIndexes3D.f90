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
 Subroutine DefineBoundaryLayerIndexes3D(Dim,NP,NF,BeginOfReg,NFR,NC,IDS,FaceType,X,Y,Z,Bound,NConnectedPoints,IConnectedPoints,NodeMetric,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,I,NP,NF,NC,NSeedNodes,LastLayerIndexNumber
 Integer,Dimension(1:Dim)::LayerIndex,ISeedNodes,FaceType,NConnectedPoints,NConnectedFaces
 Logical,Dimension(1:Dim)::StretchNode,Bound
 Real(8),Dimension(1:3,1:3,1:Dim)::NodeMetric
 Real(8),Dimension(1:Dim)::ElipsoidRatio,X,Y,Z
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints,IConnectedFaces
 Integer,Dimension(1:100)::BeginOfReg,NFR
!*********************************************************************************************
 
!Part 1:
 LayerIndex(:) = 0
 Call FindConnectedFaces(Dim,NF,NP,IDS,FaceType,IConnectedFaces,NConnectedFaces) 
 Call DetectStretchNodes3D(Dim,NP,X,Y,Z,NConnectedPoints,IConnectedPoints,NodeMetric,StretchNode,ElipsoidRatio)
 
!Part 2:
 Do I=1,NP
    If(.NOT. StretchNode(I)) LayerIndex(I)=-100000 !Layer -100000: isotropic nodes
 EndDo
 
!Part 3:
 Call DetectSeedNodes3D(Dim,NP,NF,BeginOfReg,NFR,IDS,FaceType,LayerIndex,NSeedNodes,ISeedNodes)
 
!Part 4:
 Call DoLayering3D(Dim,NF,NP,NConnectedPoints,IConnectedPoints,LayerIndex)
 
!Part 5:
 Call RemoveDiscreteLayerIndexes3D(Dim,NP,NF,IDS,FaceType,Bound,NConnectedPoints,IConnectedPoints,IConnectedFaces,NConnectedFaces,LastLayerIndexNumber,NSeedNodes,ISeedNodes,LayerIndex)
 
!Part 6:
 Call CorrectingBoundaryLayer3D(Dim,NP,NConnectedPoints,IConnectedPoints,LastLayerIndexNumber,LayerIndex)
 
!Part 7:
 Call ReDefineSeedNodes3D(Dim,NP,X,Y,Z,Bound,NConnectedPoints,IConnectedPoints,NodeMetric,LayerIndex) !Detect wake Seed Nodes
 
!Part 8:
 Call DoLayering3D(Dim,NF,NP,NConnectedPoints,IConnectedPoints,LayerIndex)
 
!Part 9:
 Call RemoveDiscreteLayerIndexes3D(Dim,NP,NF,IDS,FaceType,Bound,NConnectedPoints,IConnectedPoints,IConnectedFaces,NConnectedFaces,LastLayerIndexNumber,NSeedNodes,ISeedNodes,LayerIndex)

!*********************************************************************************************
 End
!###########################################################################################