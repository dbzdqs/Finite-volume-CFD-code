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
 Subroutine ReDefineSeedNodes3D(Dim,NP,X,Y,Z,Bound,NConnectedPoints,IConnectedPoints,NodeMetric,LayerIndex)!For Detect Wake Seed Nodes
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,Node,I,J,K,T,NewSeedNumber,MaxLayerIndex,MaxSizeNodeIndex,L1,LastNode,NCPoints,Node3,Node2,Layer
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedPoints,selectedLayerSt,ICPoints
 Real(8),Dimension(1:3,1:3,1:Dim)::NodeMetric
 Real(8)::MaxSizeValue,hSize,minVal,maxVal,maxSizeVal
 Logical::changed,changed2,firstLevel,exists,condition,firstCheck
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints
 Integer,PARAMETER:: Unknown = 0, Selected=1, OldSelected=2
 Logical,Dimension(1:Dim)::Bound
!*********************************************************************************************
!Part 1:
 selectedLayerSt(1:NP) = Unknown
 firstLevel            = .TRUE.
 changed               = .TRUE.
 
!Part 2:
 MaxLayerIndex = 0
 Do I=1,NP    
    If(LayerIndex(I)>MaxLayerIndex) MaxLayerIndex = LayerIndex(I)
 EndDo
 
!Part 3:
 Do while(changed)
    changed = .FALSE.
     
   !Find the First Seed Node And SeedIndex in each Mesh Level
   !Part 4:
    NewSeedNumber = 0
    LastNode      = 0
    Do Layer = 1,MaxLayerIndex
       Do Node = 1 , NP
          If(LayerIndex(Node)/=Layer) cycle
          If(selectedLayerSt(Node)/=Unknown) cycle
           
         !Part 5:
         !the Selected Node connected to at the selected Level
          condition = .FALSE.
          If(firstLevel .AND. Bound(Node))Then
           condition = .TRUE.
          Else If(.NOT. firstLevel)Then
           Do J=1,NConnectedPoints(Node)
              If(selectedLayerSt(IConnectedPoints(J,Node))==Selected)Then
               condition = .TRUE.
               exit
              Endif
           EndDo
          EndIf
          If(.NOT. condition) Cycle
         !the Selected Node connected to at the selected Level
            
         !Part 6:
          Do I=1,NConnectedPoints(Node)
             Node2 = IConnectedPoints(I,Node)
             If(LayerIndex(Node2)/=0)Cycle
             If(selectedLayerSt(Node2)/=Unknown) Cycle
               
            !Part 7:
            !the Selected Node be at the same Mesh Level with the first Node
             condition = .FALSE.
             If(firstLevel .AND. Bound(Node2))Then
              condition = .TRUE.
             Else If(.NOT. firstLevel)Then
              Do J=1,NConnectedPoints(Node2)
                 If(selectedLayerSt(IConnectedPoints(J,Node2))==Selected)Then
                  condition = .TRUE.
                  exit
                 Endif
              EndDo
             EndIf
             If(.NOT. condition) Cycle
            !the Selected Node be at the same Mesh Level with the first Node
                
            !Part 8:
            !the Selected Node be at the max stretching Direction From the first Node(Between LayerIndex=0 AND sameLevel Nodes)
             maxSizeVal       = 0.0
             maxSizeNodeIndex = 0
             Do J=1,NConnectedPoints(Node)
                Node3 = IConnectedPoints(J,Node)
                If(LayerIndex(Node3)/=0) Cycle
                If(selectedLayerSt(Node3)/=Unknown) Cycle
               
               !the Selected Node3 be at the same Mesh Level with the Node
                condition = .FALSE.
                If(firstLevel .AND. Bound(Node3))Then
                 condition = .TRUE.
                Else If(.NOT. firstLevel)Then
                 Do K=1,NConnectedPoints(Node3)
                    If(selectedLayerSt(IConnectedPoints(K,Node3))==Selected)Then
                     condition = .TRUE.
                     exit
                    Endif
                 EndDo
                EndIf
                If(.NOT. condition) Cycle
               !the Selected Node3 be at the same Mesh Level with the Node
               
                Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,Node),Node,Node3,hSize)
                If(hSize>maxSizeVal)Then
                 maxSizeVal       = hSize
                 maxSizeNodeIndex = Node3
                EndIf
             EndDo
             If(maxSizeNodeIndex==Node2)Then
              NewSeedNumber = Layer
              LastNode      = Node
              exit
             Endif
            !the Selected Node be at the max stretching Direction From the first Node(Between LayerIndex=0 AND sameLevel Nodes)
            
          EndDo
          If(LastNode/=0)exit
       EndDo
       If(LastNode/=0)exit
    EndDo
   !Find the First Seed Node And SeedIndex in each Mesh Level
    If(NewSeedNumber==0 .OR. LastNode==0) Cycle
      
   !Assign NewSeedNumber to new Seed Nodes in one Mesh Level
   !Part 9:
    l1 = NewSeedNumber
    changed2  = .TRUE.
    NCPoints=0
    Do while(changed2)
       changed2 = .FALSE.
         
       MaxSizeValue     = 0.0
       MaxSizeNodeIndex = 0
         
       maxSizeVal   =   0.0
       maxSizeNodeIndex = 0
       
      !Part 10:
       Do I=1,NConnectedPoints(LastNode)
          Node = IConnectedPoints(I,LastNode)
          If(LayerIndex(Node)/=0) Cycle
          If(selectedLayerSt(Node)/=Unknown) Cycle
             
         !the Selected Node3 be at the same Mesh Level with the Node
          condition = .FALSE.
          If(firstLevel .AND. Bound(Node))Then
           condition = .TRUE.
          Else If(.NOT. firstLevel)Then
           Do J=1,NConnectedPoints(Node)
              If(selectedLayerSt(IConnectedPoints(J,Node))==Selected)Then
               condition = .TRUE.
               exit
              Endif
           EndDo
          EndIf
          If(.NOT. condition) Cycle
         !the Selected Node3 be at the same Mesh Level with the Node
                  
          Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,LastNode),LastNode,Node,hSize)
          If(hSize>maxSizeVal)Then
           maxSizeVal       = hSize
           maxSizeNodeIndex = Node
          EndIf
             
       EndDo
         
       If(MaxSizeNodeIndex==0 .OR. LayerIndex(MaxSizeNodeIndex)<0) Exit
             
      !check the selected node be a stretching Node(in the current level)
      !Part 11:
       minVal     = 0.0
       maxVal     = 0.0
       firstCheck = .TRUE.
       Do J=1,NConnectedPoints(MaxSizeNodeIndex)
          Node2 = IConnectedPoints(J,MaxSizeNodeIndex)
          If(selectedLayerSt(Node2)/=Unknown) Cycle
            
          condition = .FALSE.
          If(firstLevel .AND. Bound(Node2))Then
           condition = .TRUE.
          Else If(.NOT. firstLevel)Then
           Do K=1,NConnectedPoints(Node2)
              If(selectedLayerSt(IConnectedPoints(K,Node2))==Selected)Then
               condition = .TRUE.
               exit
              Endif
           EndDo
          EndIf
          If(.NOT. condition) Cycle
            
          Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,NodeMetric(:,:,MaxSizeNodeIndex),MaxSizeNodeIndex,Node2,hSize)
          If(minVal>hSize .OR. firstCheck) minVal = hSize
          If(maxVal<hSize .OR. firstCheck) maxVal = hSize
                 
          firstCheck = .FALSE.
       EndDo
       If((maxVal/minVal)<1.60) exit
      !check the selected node be a stretching Node(in the current level)
           
      !Part 12:
       LayerIndex(MaxSizeNodeIndex) = l1
       changed  = .TRUE.
       changed2 = .TRUE.
       LastNode = MaxSizeNodeIndex
         
       exists = .FALSE.
       Do T=1,NCPoints
          If(ICPoints(T)==LastNode)Then
           exists = .TRUE.
           exit
          Endif
       EndDo
       If(.NOT. exists)Then
        NCPoints     =     NCPoints+1
        ICPoints(NCPoints) = LastNode
       Endif
         
      !Part 13:
       Do J=1,NConnectedPoints(LastNode)
          Node = IConnectedPoints(J,LastNode)
          If(selectedLayerSt(Node)/=Unknown) Cycle
            
          condition = .FALSE.
          If(firstLevel .AND. Bound(Node))Then
           condition = .TRUE.
          Else If(.NOT. firstLevel)Then
           Do K=1,NConnectedPoints(Node)
              If(selectedLayerSt(IConnectedPoints(K,Node))==Selected)Then
               condition = .TRUE.
               exit
              Endif
           EndDo
          EndIf
          If(.NOT. condition) Cycle
             
          exists = .FALSE.
          Do T=1,NCPoints
             If(ICPoints(T)==Node)Then
              exists = .TRUE.
              exit
             Endif
          EndDo
          If(.NOT. exists)Then
           NCPoints  =  NCPoints + 1
           ICPoints(NCPoints) = Node
          Endif
       EndDo  
         
    EndDo
   !Assign NewSeedNumber to new Seed Nodes in one Mesh Level
     
   !Part 14:
    If(changed)Then    
     Do I=1,NP
        If(selectedLayerSt(I)==Selected) selectedLayerSt(I) = OldSelected
     EndDo
     
    !Part 15:
     Do I=1,NCPoints
        selectedLayerSt(ICPoints(I)) = Selected
     EndDo
         
    Endif
     
    firstLevel = .FALSE.
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################