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
  Subroutine DefineBoundaryLayerIndexes(Dim,NP,NF,NC,IDS,X,Y,Bound,NConnectedEdges,IConnectedEdges,NConnectedNodes,IConnectedNodes,NodeMetric,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,I,NP,NF,NC,NSeedNodes,LastLayerIndexNumber,NExceptNodes
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedEdges,NConnectedNodes,ISeedNodes,IExceptNodes
 Logical,Dimension(1:Dim)::StretchNode,Bound,CommitedPoints
 Logical::Changed
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8),Dimension(1:Dim)::ElipsoidRatio,X,Y
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 CommitedPoints(1:NP)=.FALSE.
 Changed=.TRUE.
 NExceptNodes = 0;
 
!Part 2:
 LayerIndex(:) = 0
 Call DetectStretchNodes(Dim,NP,NodeMetric,StretchNode,ElipsoidRatio)
 
!Part 3:
 Do I=1,NP
    If(.NOT. StretchNode(I))Then
     LayerIndex(I)=-100000 !Layer -100000: isotropic nodes
     CommitedPoints(I)=.TRUE.
    Endif
 EndDo
 
 do while(Changed)
   !Part 4:
    Call DetectSeedNodes(Dim,NP,NF,IDS,Bound,ElipsoidRatio,StretchNode,NConnectedEdges,IConnectedEdges,NExceptNodes,IExceptNodes,LayerIndex,NSeedNodes,ISeedNodes,Changed)
    if(.NOT. Changed) exit
 
   !Part 5:
    Call DoLayering(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges,LayerIndex)
 
   !Part 6:
    Call RemoveDiscreteLayerIndexes(Dim,NP,NF,IDS,NConnectedEdges,IConnectedEdges,LastLayerIndexNumber,LayerIndex,CommitedPoints)
 
   !Part 7:
    Call CorrectingBoundaryLayer(Dim,NP,NConnectedNodes,IConnectedNodes,LastLayerIndexNumber,LayerIndex)
 
   !Part 8:
    Call ReDefineSeedNodes(Dim,NP,NF,IDS,X,Y,NConnectedEdges,IConnectedEdges,NodeMetric,LayerIndex) !Detect wake Seed Nodes
 
   !Part 9:
    Call DoLayering(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges,LayerIndex)
 
   !Part 10:
    Call RemoveDiscreteLayerIndexes(Dim,NP,NF,IDS,NConnectedEdges,IConnectedEdges,LastLayerIndexNumber,LayerIndex,CommitedPoints)
    
   !Part 11:
    Do I=1,NP
        If(LayerIndex(I)>0) CommitedPoints(I)=.TRUE.
    EndDo
    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################