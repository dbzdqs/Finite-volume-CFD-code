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
!// Date: Apr., 25, 2013                                                                   //!
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
 Subroutine ReDefineSeedNodes(Dim,NP,NF,IDS,X,Y,NConnectedEdges,IConnectedEdges,NodeMetric,LayerIndex)!For Detect Wake Seed Nodes
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,P2,Node,Edge,I,J,K,T,P22,NewSeedNumber,MaxLayerIndex,MaxSizeNodeIndex,L1,Edge2,LastNode
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedEdges
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric
 Real(8)::MaxSizeValue,hSize
 Logical::changed
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y
!*********************************************************************************************
!Part 1:
 MaxLayerIndex = 0
 Do I=1,NP
    If(LayerIndex(I)>MaxLayerIndex) MaxLayerIndex = LayerIndex(I)
 EndDo
 
!Part 2:
 NewSeedNumber = 0
 LastNode = 0
 Do I=1,MaxLayerIndex
    Do J=1,NP
       Node = J
       If(LayerIndex(J)/=I) Cycle
         
      !Part 3:
       Do K=1,NConnectedEdges(Node)
          Edge = IConnectedEdges(Node,K)
          If(IDS(3,Edge)==Node)Then
           P2 = IDS(4,Edge)
          Else
           P2 = IDS(3,Edge)
          Endif
          If(LayerIndex(P2)/=0)Cycle
          
         !Part 4:
          MaxSizeValue = 0.0
          MaxSizeNodeIndex = 0
          Do T=1,NConnectedEdges(Node)
             Edge2 = IConnectedEdges(Node,T)
             If(IDS(3,Edge2)==Node)Then
              P22 = IDS(4,Edge2)
             Else
              P22 = IDS(3,Edge2)
             Endif
             If(LayerIndex(P22)/=0) Cycle
                 
             Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,Node),Node,P22,hSize)
             If(hSize>MaxSizeValue)Then
              MaxSizeValue = hSize
              MaxSizeNodeIndex = P22
             Endif
          EndDo
           
         !Part 5:
          If(MaxSizeNodeIndex==P2)Then
           NewSeedNumber = LayerIndex(Node)
           LastNode = Node
           goto 10
          Endif
             
       EndDo
         
    EndDo
     
 EndDo
 
10 continue
   
!Part 6:
 If(NewSeedNumber==0 .OR. LastNode==0) return
   
!reDefine New Seed Nodes(Wake) Step By Step
!Part 7:
 l1 = NewSeedNumber
 changed = .TRUE.
 Do while (changed)
    changed = .FALSE.
    
   !Part 8:
    MaxSizeValue     = 0.0
    MaxSizeNodeIndex = 0
    Do T=1,NConnectedEdges(LastNode)
       Edge2 = IConnectedEdges(LastNode,T)
       If(IDS(3,Edge2)==LastNode)Then
        P22 = IDS(4,Edge2)
       Else
        P22 = IDS(3,Edge2)
       Endif
       If(LayerIndex(P22)>0) Cycle
           
       Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,NodeMetric(:,:,LastNode),LastNode,P22,hSize)
       If(hSize>MaxSizeValue)Then
        MaxSizeValue = hSize
        MaxSizeNodeIndex = P22
       Endif
    EndDo
    
   !Part 9:
    If(MaxSizeNodeIndex==0)exit
    If(LayerIndex(MaxSizeNodeIndex)<0)exit
    
   !Part 10:
    LayerIndex(MaxSizeNodeIndex) = l1
    changed = .TRUE.
    LastNode = MaxSizeNodeIndex
    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################