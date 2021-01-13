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
 Subroutine RemoveDiscreteLayerIndexes(Dim,NP,NF,IDS,NConnectedEdges,IConnectedEdges,LastLayerIndexNumber,LayerIndex,CommitedPoints)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,K,P22,P1,LI1,LI2,P2,Node,Edge,MaxLayerIndex,LastLayerIndexNumber
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedEdges
 Logical::exists,connectivity
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:4,1:Dim)::IDS
 Logical,Dimension(1:Dim)::CommitedPoints
!*********************************************************************************************
!Part 1:
 Do I=1,NP
    Node = I
    If(LayerIndex(Node)<=0 .OR. CommitedPoints(Node)) Cycle
     
   !Part 2:
    exists = .FALSE.
    Do J=1,NConnectedEdges(Node)
       Edge = IConnectedEdges(Node,J)
       If(IDS(3,Edge)==Node)Then
        P2 = IDS(4,Edge)
       Else
        P2 = IDS(3,Edge)
       EndIf
       If(LayerIndex(P2)>0 .AND. LayerIndex(P2)==LayerIndex(Node))Then
        exists = .TRUE.
        exit
       Endif
         
    EndDo
     
    If(.NOT. exists) LayerIndex(Node) = 0    
 EndDo
 
 
!!!! Remove Wrong Layers !!!!
!Part 3:
 MaxLayerIndex = 0
 Do I=1,NP
    If(LayerIndex(I)>MaxLayerIndex .AND. (.NOT. CommitedPoints(I))) MaxLayerIndex = LayerIndex(I)
 EndDo

!Part 4:
 LastLayerIndexNumber = 0
 Do I=1,MaxLayerIndex
    LI1 = I
    LI2 = I + 1
     
   !!!! Check Connectivity Between at least one node of the selected Layers !!!!
   !Part 5:
    connectivity = .FALSE.
    Do J=1,NP
       If(LayerIndex(J)/=LI2 .OR. CommitedPoints(J)) Cycle
       P1 = J
         
       connectivity = .FALSE.
       Do K=1,NConnectedEdges(P1)
          Edge = IConnectedEdges(P1,K)
          If(IDS(3,Edge)==P1) Then
           P22 = IDS(4,Edge)
          Else
           P22 = IDS(3,Edge)
          Endif
                 
          If(LayerIndex(P22)==LI1 .AND. (.NOT. CommitedPoints(P22)))Then
           connectivity = .TRUE.
           exit
          Endif
                 
       EndDo
       If(.NOT. connectivity) exit
    EndDo
   !!!! Check Connectivity Between at least one node of the selected Layers !!!!
   
    If(.NOT. connectivity)Then
     LastLayerIndexNumber = LI1
     exit
    Endif
     
 EndDo
 
!Part 6:
 If(LastLayerIndexNumber/=0)Then
  Do I=1,NP
     If(LayerIndex(I)>LastLayerIndexNumber .AND. (.NOT. CommitedPoints(I))) LayerIndex(I) = 0
  EndDo
 Endif
!!!! Remove Wrong Layers !!!!
!*********************************************************************************************
 End
!###########################################################################################