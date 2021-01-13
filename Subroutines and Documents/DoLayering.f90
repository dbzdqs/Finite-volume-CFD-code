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
Subroutine DoLayering(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,l1,P2,Node,Edge
 Integer,Dimension(1:Dim)::LayerIndex,NConnectedEdges
 Integer,Dimension(1:4,1:Dim)::IDS
 Logical::exists,changed
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
!*********************************************************************************************
!Part 1:
 l1 = 1
 changed = .TRUE.
 Do while (changed)
    
   !Part 2:
    changed = .FALSE.
    l1 = l1 + 1
    Do I=1,NP
       Node = I
       If(LayerIndex(Node)/=0 .OR. NConnectedEdges(Node)==0)Cycle
         
      !Part 3:
       exists = .FALSE.
       Do J=1,NConnectedEdges(Node)
          Edge = IConnectedEdges(Node,J)
           
          If(IDS(3,Edge)==Node)Then
           P2 = IDS(4,Edge)
          Else
           P2 = IDS(3,Edge)
          Endif
            
         !Part 4:
          If(LayerIndex(P2)>0 .AND. LayerIndex(P2)<l1)Then
           exists = .TRUE.
           exit
          Endif
           
       EndDo
         
      !Part 5:
       If(exists)Then
        LayerIndex(Node) = l1
        changed = .TRUE.
       Endif
       
    EndDo   
 EndDo
!*********************************************************************************************
 End
!###########################################################################################