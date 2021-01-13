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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine DetectSeedNodes(Dim,NP,NF,IDS,Bound,ElipsoidRatio,StretchNode,NConnectedEdges,IConnectedEdges,NExceptNodes,IExceptNodes,LayerIndex,NSeedNodes,ISeedNodes,hasNewSeed)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,MaxStrechNode,NSeedNodes,P2,Node,Edge,I,J,K,NExceptNodes
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedEdges,ISeedNodes,IExceptNodes
 Logical,Dimension(1:Dim)::StretchNode,Bound
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8),Dimension(1:Dim)::ElipsoidRatio
 Real(8)::MaxStrechValue
 Logical::changed,exists,hasNewSeed
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:4,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 MaxStrechValue = 0.0
 MaxStrechNode  = 0
 hasNewSeed= .FALSE.
 Do I=1,NP
    If(Bound(I) .AND. ElipsoidRatio(I)>MaxStrechValue .AND. LayerIndex(I)==0)Then
       
    !Part 2:
     exists = .FALSE.
     Do J=1,NExceptNodes
        If(IExceptNodes(J)==I)Then
         exists = .TRUE.
         exit
        Endif
     EndDo
     If(exists) Cycle
    
     MaxStrechValue = ElipsoidRatio(I)
     MaxStrechNode  = I
     hasNewSeed = .TRUE.
    Endif
 EndDo
 If(.NOT. hasNewSeed) return
 
!Part 3:
 NSeedNodes = 1
 ISeedNodes(NSeedNodes)    = MaxStrechNode
 LayerIndex(MaxStrechNode) = NSeedNodes
 NExceptNodes = NExceptNodes + 1
 IExceptNodes(NExceptNodes) = MaxStrechNode
 changed = .TRUE.
 Do while (changed)
   
   !Part 4:
    changed = .FALSE.
    Do I=1,NSeedNodes
       Node = ISeedNodes(I)
       
       Do J=1,NConnectedEdges(Node)
          Edge = IConnectedEdges(Node,J)
          If(IDS(3,Edge)==Node)Then
           P2 = IDS(4,Edge)
          Else
           P2 = IDS(3,Edge)
          Endif
           
         !Part 5:
          If(NConnectedEdges(P2)==0 .OR. (.NOT. Bound(P2)) .OR. (.NOT. StretchNode(P2)) )Cycle    
          exists=.FALSE.
          Do K=1,NSeedNodes
             If(ISeedNodes(K)==P2)Then
              exists = .TRUE.
              exit
             Endif
          EndDo
          If(.NOT. exists)Then
           NSeedNodes = NSeedNodes + 1
           ISeedNodes(NSeedNodes) = P2
           LayerIndex(P2) = 1
           NExceptNodes = NExceptNodes + 1
           IExceptNodes(NExceptNodes) = P2
           changed = .TRUE.
          Endif
           
       EndDo
    EndDo   
 EndDo
!*********************************************************************************************
 End
!###########################################################################################