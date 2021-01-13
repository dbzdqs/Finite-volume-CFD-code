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
 Subroutine RemoveDiscreteLayerIndexes3D(Dim,NP,NF,IDS,FaceType,Bound,NConnectedPoints,IConnectedPoints,IConnectedFaces,NConnectedFaces,LastLayerIndexNumber,NSeedNodes,ISeedNodes,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,K,P22,P1,LI1,LI2,P2,Node,MaxLayerIndex,LastLayerIndexNumber,Point,NCPoints,T,Point2,Face,tempP1,tempP2,NSeedNodes
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedPoints,selectedLayerSt,ICPoints,NConnectedFaces,ISeedNodes
 Logical::exists,connectivity,firstLevel,condition,changed,isConnected
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints
 Integer,PARAMETER:: Unknown = 0, Selected=1, OldSelected=2
 Integer,Dimension(1:4)::FacePoints
 Integer,Dimension(1:100,1:Dim)::IConnectedFaces
 Logical,Dimension(1:Dim)::Bound
!*********************************************************************************************
!Part 1:
 selectedLayerSt(1:NP) = Unknown
 Do I=1,NSeedNodes
     selectedLayerST(ISeedNodes(I))=OldSelected
 EndDo
 
!Part 2:
 firstLevel = .TRUE.
 changed    = .TRUE.
 Do while (changed)
    changed = .FALSE.
    Do I=1,NP   
       If(selectedLayerSt(I)/=Unknown) Cycle
        
      !Part 3:
       condition = .FALSE.
       If(firstLevel .AND. Bound(I))Then
        condition = .TRUE.
       Else If(.NOT. firstLevel)Then
        Do J=1,NConnectedPoints(I)
           If(selectedLayerSt(IConnectedPoints(J,I))==Selected)Then
            condition = .TRUE.
            exit
           Endif
        EndDo
       EndIf
       If(.NOT. condition) Cycle
        
      !Part 4:
       NCPoints = 0
       If(firstLevel)Then
        Do J=1,NConnectedFaces(I)
           Face = IConnectedFaces(J,I)
           If(IDS(2,Face)==0)Then
            FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
            Do K=1,FaceType(Face)    
               P1    = FacePoints(K)
               P2    = FacePoints(Mod(K,FaceType(Face)) + 1)
               Point = 0
               If(P1==I) Point = P2
               If(P2==I) Point = P1
               If(Point==0) Cycle
                      
               exists = .FALSE.
               Do T=1,NCPoints
                  If(ICPoints(T)==Point)Then
                   exists = .TRUE.
                   exit
                  Endif
               EndDo
                      
               If(.NOT. exists)Then
                NCPoints   =  NCPoints + 1
                ICPoints(NCPoints) = Point
               Endif
                      
            EndDo
           Endif
        EndDo
       Else
            
        Do J=1,NConnectedPoints(I)
           Point = IConnectedPoints(J,I)
                
           exists = .FALSE.
           Do T=1,NConnectedPoints(Point)
              Point2 = IConnectedPoints(T,Point)
              If(selectedLayerSt(Point2)==Selected)Then
               exists = .TRUE.
               exit
              Endif
           EndDo
           If(exists)Then
            NCPoints   =  NCPoints + 1
            ICPoints(NCPoints) = Point
           Endif
            
        EndDo
            
       Endif
           
      !Part 5:
       exists = .FALSE.
       Do J=1,NCPoints
          P2 = ICPoints(J)
          If(LayerIndex(P2)>0 .AND. LayerIndex(P2)==LayerIndex(I) .AND. selectedLayerSt(P2)==Unknown .AND. selectedLayerSt(I)==Unknown)Then  
            
           If(Bound(I) .AND. Bound(P2))Then
            exists=.TRUE.
            exit
           Endif
              
          !Part 6:
           Do K=1,NConnectedPoints(I)
              tempP1 = IConnectedPoints(K,I)
              If(selectedLayerSt(tempP1)/=Selected) Cycle
                
              Do T=1,NConnectedPoints(P2)
                 tempP2 = IConnectedPoints(T,P2)
                 If(selectedLayerSt(tempP2)/=Selected) Cycle
                    
                 isConnected=.FALSE.
                 Call CheckNodesConnectivity(tempP1,tempP2,NConnectedPoints(tempP1),IConnectedPoints(:,tempP1),NConnectedPoints(tempP2),IConnectedPoints(:,tempP2),isConnected)
                 if(isConnected)Then
                  exists = .TRUE.
                  exit
                 Endif
                    
              EndDo
              If(exists) exit
           
           EndDo       
           If(exists) exit
           
          Endif
       EndDo
       If(.NOT. exists)Then
        LayerIndex(I) = 0
        changed =  .TRUE.
       Endif
    EndDo
      
   !Part 7:
    If(changed)Then
     NCPoints=0
     Do I=1,NP
        If(firstLevel)Then
        !add all boundary points
         If(Bound(I) .AND. selectedLayerSt(I)==Unknown)Then
          NCPoints = NCPoints + 1
          ICPoints(NCPoints)  = I
         Endif
        Else    
         If(selectedLayerSt(I)/=Selected) Cycle
         Do J=1,NConnectedPoints(I)
            P2 = IConnectedPoints(J,I)
            If(selectedLayerSt(P2)==Unknown)Then
             exists = .FALSE.
             Do T=1,NCPoints
                If(ICPoints(T)==P2)Then
                 exists = .TRUE.
                 exit
                Endif
             EndDo
             If(.NOT. exists)Then
              NCPoints = NCPoints + 1
              ICPoints(NCPoints) = P2
             Endif
            Endif
         EndDo
        EndIf
             
     EndDo
         
    !Part 8:
     Do I=1,NP
        If(selectedLayerSt(I)==Selected) selectedLayerSt(I) = OldSelected
     EndDo
     Do I=1,NCPoints
        selectedLayerSt(ICPoints(I)) = Selected
     EndDo
         
    Endif
     
    firstLevel = .FALSE.    
 EndDo
 
!!!! Remove Wrong Layers !!!!
!Part 9:
 MaxLayerIndex = 0
 Do I=1,NP
    If(LayerIndex(I)>MaxLayerIndex) MaxLayerIndex = LayerIndex(I)
 EndDo

!Part 10:
 LastLayerIndexNumber = 0
 Do I=1,MaxLayerIndex
    LI1 = I
    LI2 = I + 1
     
   !!!! Check Connectivity Between at least one node of the selected Layers !!!!
   !Part 11:
    connectivity = .FALSE.
    Do J=1,NP
       If(LayerIndex(J)/=LI2) Cycle
       P1 = J
         
       connectivity = .FALSE.
       Do K=1,NConnectedPoints(P1)
          P22 = IConnectedPoints(K,P1)
              
          If(LayerIndex(P22)==LI1)Then
           connectivity = .TRUE.
           exit
          Endif
                 
       EndDo
       If(.NOT. connectivity) exit
    EndDo
   !!!! Check Connectivity Between at least one node of the selected Layers !!!!
   
   !Part 12:
    If(.NOT. connectivity)Then
     LastLayerIndexNumber = LI1
     exit
    Endif
     
 EndDo
 
!Part 13:
 If(LastLayerIndexNumber/=0)Then
  Do I=1,NP
     If(LayerIndex(I)>LastLayerIndexNumber) LayerIndex(I) = 0
  EndDo
 Endif
!!!! Remove Wrong Layers !!!!
!*********************************************************************************************
 End
!###########################################################################################