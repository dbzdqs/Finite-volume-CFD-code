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
 Subroutine DoLayering3D(Dim,NF,NP,NConnectedPoints,IConnectedPoints,LayerIndex)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,l1,P2,Node
 Integer,Dimension(1:Dim)::LayerIndex
 Logical::exists,changed
 Integer,Dimension(1:Dim)::NConnectedPoints
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints
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
       If(LayerIndex(Node)/=0 .OR. NConnectedPoints(Node)==0)Cycle
         
      !Part 3:
       exists = .FALSE.
       Do J=1,NConnectedPoints(Node)
          P2 = IConnectedPoints(J,Node)
           
            
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