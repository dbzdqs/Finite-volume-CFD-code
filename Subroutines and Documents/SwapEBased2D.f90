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
!// Date: Oct., 05, 2016                                                                   //!
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
 Subroutine SwapEBased2D(Dim,EI,E1,E2,E3,E4,P1,P2,P3,P4,IDS,InxEdgeOfCell,Stack,NStack)
 Implicit None
!*********************************************************************************************
 Intent(In     )::Dim,EI,E1,E2,E3,E4,P1,P2,P3,P4 
 Intent(Inout  )::InxEdgeOfCell,Stack,NStack
 
 Integer::Dim,ME,NE,EI,E1,E2,E3,E4,P1,P2,P3,P4,NStack
 Integer,Dimension(1:4, 1:Dim)::IDS
 Integer,Dimension(1:4, 1:Dim)::InxEdgeOfCell
 Integer,Dimension(1:Dim)::Stack
!*********************************************************************************************
!Part 1:
 ME = IDS(1, EI)
 NE = IDS(2, EI)

!Part 2:
 IDS(3, EI) = P3
 IDS(4, EI) = P4

!Part 3:
 If(IDS(1, E1) == ME) Then
  IDS(1, E1) = NE
 Elseif(IDS(2, E1) == ME) Then
  IDS(2, E1) = NE
 Endif

!Part 4:
 If(IDS(1, E2) == NE) Then
  IDS(1, E2) = ME
 Elseif(IDS(2, E2) == NE) Then
  IDS(2, E2) = ME
 Endif 

!Part 5:    
 InxEdgeOfCell(1,ME) = E3
 InxEdgeOfCell(2,ME) = EI
 InxEdgeOfCell(3,ME) = E2

!Part 6:  
 InxEdgeOfCell(1,NE) = E4
 InxEdgeOfCell(2,NE) = EI
 InxEdgeOfCell(3,NE) = E1  

!Part 7:
 If( IDS(1, E1)/=0 .and.  IDS(2, E1)/=0 ) Then
  NStack = NStack + 1
  Stack(NStack) = E1
 Endif

 If( IDS(1, E2)/=0 .and.  IDS(2, E2)/=0 ) Then
  NStack = NStack + 1
  Stack(NStack) = E2
 Endif

 If( IDS(1, E3)/=0 .and.  IDS(2, E3)/=0 ) Then
  NStack = NStack + 1
  Stack(NStack) = E3
 Endif

 If( IDS(1, E4)/=0 .and.  IDS(2, E4)/=0 ) Then
  NStack = NStack + 1
  Stack(NStack) = E4
 Endif
!*********************************************************************************************
 End
!###########################################################################################