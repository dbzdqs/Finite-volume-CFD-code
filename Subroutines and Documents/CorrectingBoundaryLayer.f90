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
!// Date: Nov., 15, 2014                                                                   //!
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
 Subroutine CorrectingBoundaryLayer(Dim,NP,NConnectedNodes,IConnectedNodes,LastLayerIndexNumber,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,L,P2,Node,LastLayerIndexNumber,Layer
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedNodes
 Logical::exists1,exists2,changed
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
!*********************************************************************************************
!Part 1:
 Do L=1,LastLayerIndexNumber
    Layer = L
    
   !Part 2:
    changed = .TRUE.
    Do while (changed)
    
      !Part 3:
       changed = .FALSE.
       Do I=1,NP
          Node = I 
          If(LayerIndex(Node)>0 .OR. NConnectedNodes(Node)==0) Cycle
         
         !Part 4:
          exists1 = .FALSE.
          exists2 = .FALSE.
          Do J=1,NConnectedNodes(Node)
             P2 = IConnectedNodes(J,Node)
             If(abs(LayerIndex(P2))>0 .AND. (abs(LayerIndex(P2))==(abs(Layer)-1))  ) exists1 = .TRUE.
             If(abs(LayerIndex(P2))>0 .AND. (abs(LayerIndex(P2))==(abs(Layer)))) exists2 = .TRUE.
          EndDo
         
         !Part 5:
          If(exists1 .AND. exists2)Then
           LayerIndex(Node) = Layer
           changed = .TRUE.
          Endif
       
       EndDo   
    EndDo    
 EndDo
!*********************************************************************************************
 End
!###########################################################################################