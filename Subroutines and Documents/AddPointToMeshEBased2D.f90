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
 Subroutine AddPointToMeshEBased2D(Dim,IDS,NP,NF,NC,CandCell,X,Y,InxEdgeOfCell,NStack,Stack)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,CandCell
 Intent(InOut)::IDS,NF,NC,X,Y,InxEdgeOfCell,NStack,Stack

 Integer::Dim,NP,NC,NF,CandCell,NStack,E1,E2,E3,P1,P2,P3
 Integer,Dimension(1:4, 1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X, Y
 Integer,Dimension(1:4,1:Dim)::InxEdgeOfCell
 Integer,Dimension(1:Dim)::Stack
!*********************************************************************************************
!Part 1:
 E1 = InxEdgeOfCell(1, CandCell  )
 E2 = InxEdgeOfCell(2, CandCell  )
 E3 = InxEdgeOfCell(3, CandCell  )

!Part 2:
 IF( E1>0 )Then
  P1 = IDS(3,E1)
  P2 = IDS(4,E1)
 Else
  P1 = IDS(4,-E1)
  P2 = IDS(3,-E1)
 EndIF

!Part 3:
 IF( E2>0 )Then
  P3 = IDS(4,E2)
 Else
  P3 = IDS(3,-E2)
 EndIF

!Part 4:
 IDS(1, NF+1) = NC+2
 IDS(2, NF+1) = CandCell  
 IDS(3, NF+1) = P1
 IDS(4, NF+1) = NP

!Part 5:
 IDS(1, NF+2) = CandCell  
 IDS(2, NF+2) = NC+1
 IDS(3, NF+2) = P2
 IDS(4, NF+2) = NP

!Part 6:
 IDS(1, NF+3) = NC+1
 IDS(2, NF+3) = NC+2
 IDS(3, NF+3) = P3
 IDS(4, NF+3) = NP

!Part 7:
 IF( E2>0 )Then
  IDS(1, E2)  = NC+1
 Else
  IDS(2, -E2) = NC+1
 EndIF

!Part 8:
 IF( E3>0 )Then
  IDS(1, E3)  = NC+2
 Else
  IDS(2, -E3) = NC+2
 EndIF

!Part 9:
 InxEdgeOfCell(1,CandCell  ) = E1
 InxEdgeOfCell(2,CandCell  ) = NF+2
 InxEdgeOfCell(3,CandCell  ) =-(NF+1)

!Part 10:
 InxEdgeOfCell(1,NC+1) = E2
 InxEdgeOfCell(2,NC+1) = NF+3
 InxEdgeOfCell(3,NC+1) =-(NF+2)

!Part 11:
 InxEdgeOfCell(1,NC+2) = E3
 InxEdgeOfCell(2,NC+2) = NF+1
 InxEdgeOfCell(3,NC+2) =-(NF+3)

!Part 12:
 NStack = NStack + 1 
 Stack(NStack) = NF+1
 
 NStack = NStack + 1
 Stack(NStack) = NF+2

 NStack = NStack + 1
 Stack(NStack) = NF+3 

 NStack = NStack + 1 
 Stack(NStack) = E1
 
 NStack = NStack + 1
 Stack(NStack) = E2
  
 NStack = NStack + 1
 Stack(NStack) = E3

!Part 13:
 NC = NC + 2
 NF = NF + 3
!*********************************************************************************************
 End
!###########################################################################################