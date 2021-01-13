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
Subroutine CheckNeighborVertices(Dim,NP,NF,NC,Bound,Vertice,NConnectedNodes,IConnectedNodes,VNN,NodeMetric,CoarsedNodeMetric,X,Y,LayerIndex,VST,NContractedEdges,IContractedEdges)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,I,J,T,J2,K,Vertice,NeibV,NeibV2,NContractedEdges,NConnectedVertices,counter,Node,NodeIndex,NCheckedNodes
 Real(8),Dimension(1:Dim)::X,Y,VNN
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
 Logical::exists,ControlFaild
 Real(8),Dimension(1:2,1:2,1:Dim)::NodeMetric,CoarsedNodeMetric
 Real(8)::hV,hNV,EdgeLength,EdgeLength2,threshold,tmpX,tmpY
 Logical,Dimension(1:Dim)::Bound
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedNodes,VST !Vertice Status: 0:unknown - 1:included - 2:excluded
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
 Integer,PARAMETER:: Unknown = 0, Included=1, Excluded=2
 Integer,Dimension(200)::ICheckedNodes
!*********************************************************************************************

!Part 1:
 NConnectedVertices=NConnectedNodes(Vertice)
 counter=1
 NCheckedNodes=0
 Do while (counter<=NConnectedVertices)
    counter = counter + 1
    
   !Part 2:
    Call FindMinStretchNeib(Dim,NP,Vertice,NConnectedNodes(Vertice),IConnectedNodes(:,Vertice),NodeMetric,NCheckedNodes,ICheckedNodes,NodeIndex)
    If(NodeIndex==0) Cycle
    NeibV = NodeIndex
     
    NCheckedNodes = NCheckedNodes + 1
    ICheckedNodes(NCheckedNodes) = NeibV
     
   !Part 3:
    If(VST(NeibV)==Excluded)Cycle    
    If(LayerIndex(NeibV)>0 .AND. LayerIndex(Vertice)<=0) Cycle
    If(Bound(NeibV) .AND. ( .NOT. (Bound(Vertice))))Cycle
     
   !Part 4:
    EdgeLength = sqrt(( (X(NeibV)-X(Vertice)) **2)+( (Y(NeibV)-Y(Vertice)) **2))
    If(LayerIndex(NeibV)<=0)Then
     If( (sqrt(2.0)*(VNN(Vertice)+VNN(NeibV))) < ((EdgeLength)) ) Cycle
    Else 
     Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,CoarsedNodeMetric(:,:,NeibV),NeibV,Vertice,hNV)
     If( ((hNV)) < ((EdgeLength)) ) Cycle    
    Endif
     
   !!! Control Mesh Gradation !!!
   !Part 5:
    tmpX=X(NeibV)
    tmpY=Y(NeibV)
    X(NeibV)=X(Vertice)
    Y(NeibV)=Y(Vertice)
    ControlFaild = .FALSE.
    Do J=1,NConnectedNodes(NeibV)
       Node = IConnectedNodes(J,NeibV)
       If(Node==Vertice) Cycle     
       EdgeLength2=sqrt(( (X(Node)-X(Vertice)) **2)+( (Y(Node)-Y(Vertice)) **2))
            
       If(LayerIndex(NeibV)<=0)Then
        threshold = (sqrt(2.0)*abs(VNN(Vertice)+VNN(Node)))/(EdgeLength2)
       Else
        Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,CoarsedNodeMetric(:,:,Node),Node,Vertice,hNV)
        Call CalcLocalSizeInNodeDirection(Dim,NP,X,Y,CoarsedNodeMetric(:,:,Vertice),Vertice,Node,hV)
        threshold = ((hNV+hV)/EdgeLength2)
       Endif
         
       If(threshold<=0.5)Then
        ControlFaild = .TRUE.
        exit
       Endif
            
    EndDo
    X(NeibV)=tmpX
    Y(NeibV)=tmpY
   !!! Control Mesh Gradation !!!
    If(ControlFaild) cycle
     
   !Part 6:
    VST(Vertice) = Included
    VST(NeibV)   = Excluded
     
    NContractedEdges = NContractedEdges + 1
    IContractedEdges(1,NContractedEdges) = Vertice   !Heir
    IContractedEdges(2,NContractedEdges) = NeibV !Dead
     
   !Part 7:
   !!!! Remove ExNode From Vertice Neibs !!!!
    Do J=1,NConnectedNodes(Vertice)        
       If(IConnectedNodes(J,Vertice)==NeibV .AND. J<=NConnectedNodes(Vertice))Then
        IConnectedNodes(J,Vertice) = IConnectedNodes(NConnectedNodes(Vertice),Vertice)
        NConnectedNodes(Vertice) = NConnectedNodes(Vertice) - 1
        counter = counter - 1
        exit
       Endif
    EndDo
   !!!! Remove ExNode From Vertice Neibs !!!!
     
   !!!! Add ExNode Neibs To Vertice Neibs AND Vertice To ExNode Neibs !!!!
   !Part 8:
    Do J=1,NConnectedNodes(NeibV)
       NeibV2 = IConnectedNodes(J,NEibV)
       If(NeibV2==Vertice) Cycle
         
       exists = .FALSE.
       Do K=1,NConnectedNodes(Vertice)
          If(IConnectedNodes(K,Vertice)==NeibV2)Then
           exists = .TRUE.
           exit
          Endif
       EndDo
       If(.NOT. exists)Then
        NConnectedNodes(Vertice) = NConnectedNodes(Vertice) + 1
        IConnectedNodes(NConnectedNodes(Vertice),Vertice) = NeibV2
       Endif
         
       Do K=1,NConnectedNodes(NeibV2)
             
          If(IConnectedNodes(K,NeibV2)==NeibV)Then
           exists = .FALSE.
           Do T=1,NConnectedNodes(NeibV2)
              If(IConnectedNodes(T,NeibV2)==Vertice)Then
               exists = .TRUE.
               exit
              Endif
           EndDo
           If(.NOT. exists)Then
            IConnectedNodes(K,NeibV2) = Vertice
           Else
            IConnectedNodes(K,NeibV2) = IConnectedNodes(NConnectedNodes(NeibV2),NeibV2)
            NConnectedNodes(NeibV2) = NConnectedNodes(NeibV2) - 1
            exit
           Endif
          Endif
       EndDo
         
    EndDo 
   !!!! Add ExNode Neibs To Vertice Neibs AND Vertice To ExNode Neibs !!!!
     
   !Part 9:
    NConnectedVertices = NConnectedNodes(Vertice)
    If(LayerIndex(Vertice)>1) exit
     
 EndDo
 
!*********************************************************************************************
 End
!###########################################################################################