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
 Subroutine CheckNeighborVertices3D(Dim,NP,NF,NC,IDS,FaceType,NConectCell,IConectCell,NFace_Cell,IFace_Cell,Bound,Vertice,NConnectedPoints,IConnectedPoints,VNN,NodeMetric,CoarsedNodeMetric,X,Y,Z,LayerIndex,VST,NContractedEdges,IContractedEdges)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NC,I,J,T,J2,K,Vertice,NeibV,NeibV2,NContractedEdges,NConnectedVertices,counter,Node,NodeIndex,NCheckedNodes,MinNodeIndex,MaxNodeIndex,AllowedLevelNum,minDirection,maxDirection2,ierr
 Real(8)::hV,hNV,EdgeLength,EdgeLength2,threshold,tmpX,tmpY,tmpZ,MinValue,MaxValue,minValueDir,maxValue2
 Integer,Dimension(1:Dim)::NConnectedPoints,FaceType,NFace_Cell,NConectCell,LayerIndex,VST !Vertice Status: 0:unknown - 1:included - 2:excluded
 Real(8),Dimension(1:3,1:3,1:Dim)::NodeMetric,CoarsedNodeMetric
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints,IConectCell
 Integer,PARAMETER:: Unknown = 0, Included=1, Excluded=2
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
 Real(8),Dimension(1:Dim)::X,Y,Z,VNN
 Real(8),Dimension(1:3)::EigenValue,localSize
 Logical::exists,ControlFaild
 Logical,Dimension(1:Dim)::Bound
 Integer,Dimension(100)::ICheckedNodes
 Real(8),Dimension(3)::coordinateDiff
 Real(8),Dimension(1:3,1:3)::U,Vt
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
!*********************************************************************************************
!Part 1:
 CALL SVD(NodeMetric(:,:,Vertice),3, EigenValue, .TRUE., U, .TRUE., Vt, ierr)
 localSize(:)=sqrt(1.0/EigenValue(:))
 minDirection=0
 minValueDir=0.0
 Do J=1,3
     If(minValueDir>localSize(J) .OR. minDirection==0)Then
         minDirection=J
         minValueDir=localSize(J)
     Endif
 EndDo
 
!Part 2:
 NConnectedVertices = NConnectedPoints(Vertice)
 counter=1
 NCheckedNodes=0
 Do while (counter<=NConnectedVertices)
    counter = counter + 1
    
   !Part 3:
    Call FindMinMaxStretchNeib(Dim,NP,Vertice,X,Y,Z,NConnectedPoints(Vertice),IConnectedPoints(:,Vertice),NodeMetric,NCheckedNodes,ICheckedNodes,MinNodeIndex,MaxNodeIndex,MinValue,MaxValue)
    If(MinNodeIndex==0) Cycle
    NeibV = MinNodeIndex
    
    NCheckedNodes = NCheckedNodes + 1
    ICheckedNodes(NCheckedNodes) = NeibV
     
   !Part 4:
    If(VST(NeibV)==Excluded .OR. VST(NeibV)==Included) cycle
    If(LayerIndex(NeibV)>0 .AND. LayerIndex(Vertice)<=0) cycle
    If(Bound(NeibV) .AND. ( .NOT. (Bound(Vertice)))) cycle
    
   !Part 5:
    EdgeLength = sqrt(( (X(NeibV)-X(Vertice)) **2)+( (Y(NeibV)-Y(Vertice)) **2)+( (Z(NeibV)-Z(Vertice)) **2))
    If(LayerIndex(NeibV)<=0)Then
    !Part 6:
     If( (sqrt(2.0)*(VNN(Vertice)+VNN(NeibV))) < ((EdgeLength)) ) Cycle
    Else 
        
    !Part 7:
     Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,CoarsedNodeMetric(:,:,Vertice),Vertice,NeibV,hV)
     If(hV<EdgeLength )cycle
    
     coordinateDiff(1)=abs(X(NeibV)-X(Vertice))
     coordinateDiff(2)=abs(Y(NeibV)-Y(Vertice))
     coordinateDiff(3)=abs(Z(NeibV)-Z(Vertice))
     maxDirection2=0
     maxValue2=0.0
     Do J=1,3
         If(maxValue2<coordinateDiff(J) .OR. maxDirection2==0)Then
             maxDirection2=J
             maxValue2=coordinateDiff(J)
         Endif
     EndDo
     If(maxDirection2/=minDirection) cycle
     
    Endif
     
   !!! Control Mesh Gradation !!!
   !Part 8:
    tmpX=X(NeibV)
    tmpY=Y(NeibV)
    tmpZ=Z(NeibV)
    X(NeibV)=X(Vertice)
    Y(NeibV)=Y(Vertice)
    Z(NeibV)=Z(Vertice)
    ControlFaild = .FALSE.
    Do J=1,NConnectedPoints(NeibV)
       Node = IConnectedPoints(J,NeibV)
       If(Node==Vertice) Cycle     
       EdgeLength2=sqrt(( (X(Node)-X(Vertice)) **2)+( (Y(Node)-Y(Vertice)) **2)+( (Z(Node)-Z(Vertice)) **2))
           
      !Part 9:
       If(LayerIndex(NeibV)<=0)Then
        threshold = (sqrt(2.0)*abs(VNN(Vertice)+VNN(Node)))/(EdgeLength2)
       Else
        Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,CoarsedNodeMetric(:,:,Node),Node,Vertice,hNV)
        Call CalcLocalSizeInNodeDirection3D(Dim,NP,X,Y,Z,CoarsedNodeMetric(:,:,Vertice),Vertice,Node,hV)
        threshold = ((hNV+hV)/EdgeLength2)
       Endif
         
      !Part 10:
       If(threshold<=0.5)Then
        ControlFaild = .TRUE.
        exit
       Endif
            
    EndDo
    
   !Part 11:
    If(.NOT. ControlFaild)Then
        Call VolumeCheck3D(Dim,NF,NC,NP,IDS,X,Y,Z,FaceType,NConectCell(NeibV),IConectCell(:,NeibV),NFace_Cell,IFace_Cell,ControlFaild)
    Endif
    
   !Part 12:
    X(NeibV)=tmpX
    Y(NeibV)=tmpY
    Z(NeibV)=tmpZ
    
    If(ControlFaild) cycle
   !!! Control Mesh Gradation !!!
     
   !Part 13:
    VST(Vertice) = Included
    VST(NeibV)   = Excluded
     
    NContractedEdges = NContractedEdges + 1
    IContractedEdges(1,NContractedEdges) = Vertice  !Heir
    IContractedEdges(2,NContractedEdges) = NeibV    !Dead
     
   !Part 14:
   !!!! Remove ExNode From Vertice Neibs !!!!
    Do J=1,NConnectedPoints(Vertice)        
       If(IConnectedPoints(J,Vertice)==NeibV .AND. J<=NConnectedPoints(Vertice))Then
        IConnectedPoints(J,Vertice) = IConnectedPoints(NConnectedPoints(Vertice),Vertice)
        NConnectedPoints(Vertice) = NConnectedPoints(Vertice) - 1
        counter = counter - 1
        exit
       Endif
    EndDo
   !!!! Remove ExNode From Vertice Neibs !!!!
     
   !!!! Add ExNode Neibs To Vertice Neibs AND Vertice To ExNode Neibs !!!!
   !Part 15:
    Do J=1,NConnectedPoints(NeibV)
       NeibV2 = IConnectedPoints(J,NeibV)
       If(NeibV2==Vertice) Cycle
         
       exists = .FALSE.
       Do K=1,NConnectedPoints(Vertice)
          If(IConnectedPoints(K,Vertice)==NeibV2)Then
           exists = .TRUE.
           exit
          Endif
       EndDo
       If(.NOT. exists)Then
        NConnectedPoints(Vertice) = NConnectedPoints(Vertice) + 1
        IConnectedPoints(NConnectedPoints(Vertice),Vertice) = NeibV2
       Endif
         
      !Part 16:
       Do K=1,NConnectedPoints(NeibV2)
             
          If(IConnectedPoints(K,NeibV2)==NeibV)Then
           exists = .FALSE.
           Do T=1,NConnectedPoints(NeibV2)
              If(IConnectedPoints(T,NeibV2)==Vertice)Then
               exists = .TRUE.
               exit
              Endif
           EndDo
           If(.NOT. exists)Then
            IConnectedPoints(K,NeibV2) = Vertice
           Else
            IConnectedPoints(K,NeibV2) = IConnectedPoints(NConnectedPoints(NeibV2),NeibV2)
            NConnectedPoints(NeibV2) = NConnectedPoints(NeibV2) - 1
            exit
           Endif
          Endif
       EndDo
         
    EndDo 
   !!!! Add ExNode Neibs To Vertice Neibs AND Vertice To ExNode Neibs !!!!
     
   !Part 17:
    NConnectedVertices = NConnectedPoints(Vertice)
    If(LayerIndex(Vertice)>1) exit
     
 EndDo
 
!*********************************************************************************************
 End
!###########################################################################################