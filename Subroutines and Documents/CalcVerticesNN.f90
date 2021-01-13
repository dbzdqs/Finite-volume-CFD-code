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
 Subroutine CalcVerticesNN(Dim,NP,X,Y,NConnectedNodes,IConnectedNodes,VNN)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,I,J,Node1,Node2
 Real(8),Dimension(1:Dim)::X,Y,VNN
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
 Integer,Dimension(1:Dim)::NConnectedNodes
 Real(8)::MinDist,hV
!*********************************************************************************************
!Part 1:
 Do I=1,NP    
    Node1 = I
    If (NConnectedNodes(Node1)==0) Cycle
    MinDist  = 0
    
   !Part 2:
    Do J=1,NConnectedNodes(Node1)
       Node2 = IConnectedNodes(J,Node1)
       hV=sqrt(((X(Node2)-X(Node1))**2)+((Y(Node2)-Y(Node1))**2))
       
      !Part 5:
       If(hV < MinDist .OR. MinDist==0) MinDist = hV  
    EndDo
     
   !Part 6:
    VNN(I) = (MinDist/2.0)
 EndDo
!*********************************************************************************************
 End
!###########################################################################################