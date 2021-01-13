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
!// Developed by: F. Farhadkhani, Mathmatical, Amirkabir university of Technology          //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Do_DelaunayEBased2D(Dim,IDS,X,Y,InxEdgeOfCell,Stack,NStack)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,X,Y
 Intent(InOut)::IDS,InxEdgeOfCell,Stack,NStack

 Integer::Dim, EI,P3,P4, E1, E2, E3, E4, p1, p2, NStack,DeL
 Integer,Dimension(1:4, 1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X, Y
 Integer,Dimension(1:Dim)::Stack
 Integer,Dimension(1:4, 1:Dim)::InxEdgeOfCell
!*********************************************************************************************
!Part 1:
 Do while(NStack /= 0)

   !Part 2:
    EI = Abs( Stack(NStack) )
    NStack = NStack - 1

    IF(IDS(2,EI)==0)Cycle

   !Part 3:
    P1 = IDS(3, EI)
    P2 = IDS(4, EI)

   !Part 4:
    call Find_PointEdge(Dim,EI,IDS,InxEdgeOfCell,P3,P4,E1,E2,E3,E4)

   !Part 5:
    call DelCeckEBased2D(Dim,P1,P2,P3,P4,X,Y,DeL)

   !Part 6:
    If(DeL == -1) call SwapEBased2D(Dim,EI,E1,E2,E3,E4,P1,P2,P3,P4,IDS,InxEdgeOfCell,Stack,NStack)

 End do
!*********************************************************************************************
 End
!###########################################################################################
