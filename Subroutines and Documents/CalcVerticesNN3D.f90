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
 Subroutine CalcVerticesNN3D(Dim,NP,X,Y,Z,NConnectedPoints,IConnectedPoints,VNN)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,Node1,Node2
 Real(8),Dimension(1:Dim)::X,Y,Z,VNN
 Real(8)::MinDist,hV
 Integer,Dimension(1:Dim)::NConnectedPoints
 Integer,Dimension(1:100,1:Dim)::IConnectedPoints
!*********************************************************************************************
!Part 1:
 Do I=1,NP    
    Node1 = I
    If (NConnectedPoints(Node1)==0) Cycle
    MinDist  = 0
    
   !Part 2:
    Do J=1,NConnectedPoints(Node1)
       Node2 = IConnectedPoints(J,Node1)
       
      !Part 3:
       hV=sqrt(((X(Node2)-X(Node1))**2)+((Y(Node2)-Y(Node1))**2)+((Z(Node2)-Z(Node1))**2))
       If(hV < MinDist .OR. MinDist==0) MinDist = hV  
    EndDo
     
   !Part 4:
    VNN(I) = (MinDist/2.0)
 EndDo
!*********************************************************************************************
 End
!###########################################################################################