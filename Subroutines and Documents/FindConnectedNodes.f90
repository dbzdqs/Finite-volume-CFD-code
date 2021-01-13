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
 Subroutine FindConnectedNodes(Dim,NF,NP,IDS,NConnectedEdges,IConnectedEdges,NConnectedNodes,IConnectedNodes)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,I,J,P1,P2,Edge
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim,1:100)::IConnectedEdges
 Integer,Dimension(1:100,1:Dim)::IConnectedNodes
 Integer,Dimension(1:Dim)::NConnectedEdges,NConnectedNodes
!*********************************************************************************************
!Part 1:
 NConnectedNodes(1:NP) = 0
 Do I=1,NP
    P1 = I
     
   !Part 2:
    Do J=1,NConnectedEdges(P1)
       Edge = IConnectedEdges(P1,J)
       If(IDS(3,Edge)==P1)Then
        P2 = IDS(4,Edge)
       Else
        P2 = IDS(3,Edge)
       Endif
         
      !Part 3:
       NConnectedNodes(P1) = NConnectedNodes(P1) + 1
       IConnectedNodes(NConnectedNodes(P1),P1)   = P2
         
    EndDo
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################