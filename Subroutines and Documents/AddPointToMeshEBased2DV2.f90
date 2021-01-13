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
!// Date: May., 15, 2016                                                                   //!
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
 Subroutine AddPointToMeshEBased2DV2(Dim,IDS,NP,NF,NC,CandCell,X,Y,CELL_EDGE,NStack,Stack,NFR,InxReg)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NP,CandCell
 Intent(InOut)::IDS,NF,NC,X,Y,CELL_EDGE,NStack,Stack

 Integer::Dim,NP,NC,NF,CandCell,NStack,E1,E2,E3,P1,P2,P3,LastFacOfReg,InxReg,I,J
 Integer,Dimension(1:4, 1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X, Y
 Integer,Dimension(1:4,1:Dim)::CELL_EDGE
 Integer,Dimension(1:Dim)::Stack
 Integer,Dimension(1:100)::NFR
!*********************************************************************************************
!Part 1:
 E1 = CELL_EDGE(1, CandCell)
 E2 = CELL_EDGE(2, CandCell)
 E3 = CELL_EDGE(3, CandCell)

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
 LastFacOfReg = 0
 
!Part 5:
 Do I=1,InxReg
    LastFacOfReg = LastFacOfReg + NFR(I)
 End Do
 
!Part 6:
 Do I=NF,LastFacOfReg+1,-1
    IDS(:,I+3) = IDS(:,I)
 End Do
 
!Part 7:
 Do I=1,NC
    Do J=1,3 
       IF( CELL_EDGE(J,I)>LastFacOfReg ) CELL_EDGE(J,I)=CELL_EDGE(J,I)+3
    End Do
 End Do
 
!Part 8: 
 E1 = CELL_EDGE(1, CandCell)
 E2 = CELL_EDGE(2, CandCell)
 E3 = CELL_EDGE(3, CandCell)
 
!Part 9:
 IDS(1, LastFacOfReg+1) = NC+2
 IDS(2, LastFacOfReg+1) = CandCell
 IDS(3, LastFacOfReg+1) = P1
 IDS(4, LastFacOfReg+1) = NP

!Part 10:
 IDS(1, LastFacOfReg+2) = CandCell
 IDS(2, LastFacOfReg+2) = NC+1
 IDS(3, LastFacOfReg+2) = P2
 IDS(4, LastFacOfReg+2) = NP

!Part 11:
 IDS(1, LastFacOfReg+3) = NC+1
 IDS(2, LastFacOfReg+3) = NC+2
 IDS(3, LastFacOfReg+3) = P3
 IDS(4, LastFacOfReg+3) = NP

!Part 12:
 IF( E2>0 )Then
  IDS(1, E2)  = NC+1
 Else
  IDS(2, -E2) = NC+1
 EndIF

!Part 13:
 IF( E3>0 )Then
  IDS(1, E3)  = NC+2
 Else
  IDS(2, -E3) = NC+2
 EndIF

!Part 14:
 CELL_EDGE(1,CandCell) = E1
 CELL_EDGE(2,CandCell) = LastFacOfReg+2
 CELL_EDGE(3,CandCell) =-(LastFacOfReg+1)

!Part 15:
 CELL_EDGE(1,NC+1) = E2
 CELL_EDGE(2,NC+1) = LastFacOfReg+3
 CELL_EDGE(3,NC+1) =-(LastFacOfReg+2)

!Part 16:
 CELL_EDGE(1,NC+2) = E3
 CELL_EDGE(2,NC+2) = LastFacOfReg+1
 CELL_EDGE(3,NC+2) =-(LastFacOfReg+3)

!Part 17:
 NStack = NStack + 1 
 Stack(NStack) = LastFacOfReg+1
 
 NStack = NStack + 1
 Stack(NStack) = LastFacOfReg+2

 NStack = NStack + 1
 Stack(NStack) = LastFacOfReg+3 

 NStack = NStack + 1 
 Stack(NStack) = E1
 
 NStack = NStack + 1
 Stack(NStack) = E2
  
 NStack = NStack + 1
 Stack(NStack) = E3

!Part 18:
 NC = NC + 2
 NF = NF + 3
 NFR(InxReg) = NFR(InxReg) + 3
!*********************************************************************************************
 End
!###########################################################################################
