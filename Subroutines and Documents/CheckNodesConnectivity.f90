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
 Subroutine CheckNodesConnectivity(Node1,Node2,P1NConnectedPoints,P1IConnectedPoints,P2NConnectedPoints,P2IConnectedPoints,isConnected)
 Implicit None
!*********************************************************************************************
 Integer::Node1,Node2,I,J,P1,P2,P1NConnectedPoints,P2NConnectedPoints
 Logical::isConnected
 Integer,Dimension(1:100)::P1IConnectedPoints,P2IConnectedPoints
!*********************************************************************************************
!Part 1:
 isConnected = .FALSE.
 Do I=1,P1NConnectedPoints
    P1 = P1IConnectedPoints(I)
    
   !Part 2:
    If(P1==Node2)Then
     isConnected = .TRUE.
     exit
    Endif
     
   !Part 3:
    Do J=1,P2NConnectedPoints
       P2 = P2IConnectedPoints(J)
       If(P2==Node1)Then
        isConnected = .TRUE.
        exit
       Endif
        
      !Part 4:
       if(P1==P2)Then
        isConnected = .TRUE.
        exit
       Endif
        
    EndDo
    if(isConnected)exit
 EndDo
!*********************************************************************************************
 End
!###########################################################################################