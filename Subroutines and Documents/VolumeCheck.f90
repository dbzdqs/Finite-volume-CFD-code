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
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine VolumeCheck(Dim,Dead,heir,NConectCell,IConectCell,IEdgeOfCell,NEdgeOfCell,IDS,X,Y,NegativeVol)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Dead,heir,NConectCell,IConectCell,IEdgeOfCell,NEdgeOfCell,IDS
 Intent(InOut)::X,Y,NegativeVol

 Integer::Dim,J,JJ,E,P1,P2,NConectCell,Cell,NegativeVol,Dead,heir
 Real(8)::TmpX,TmpY,Area,DArea
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::IConectCell
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(4,DIM)::IEdgeOfCell
 Real(8),Dimension(1:Dim)::X,Y
 
 Logical::hasDeadEdge
!*********************************************************************************************
!Part 1:
 NegativeVol = -1

!Part 2:
 TmpX    = X(Dead)
 TmpY    = Y(Dead)

!Part 3:
 X(Dead) = X(heir)
 Y(Dead) = Y(heir)

!Part 4:
 Do J=1,NConectCell
    
   !Part 5:
    Cell = IConectCell(J)
    hasDeadEdge = .FALSE.

   !Part 6:
    Area = 0.0
    Do JJ=1,NEdgeOfCell(Cell)
       E = IEdgeOfCell(JJ,Cell)
       
      !Part 7:
       IF( IDS(1,E)==Cell )Then
        P1 = IDS(3,E)
        P2 = IDS(4,E)
       Else
        P1 = IDS(4,E)
        P2 = IDS(3,E)
       EndIF
       If((P1==Dead .AND. P2==Heir) .OR. (P1==Heir .AND. P2==Dead))Then
        hasDeadEdge = .TRUE.
       Endif
       
       
      !Part 8:
       DArea    = X(P1)*Y(P2) - X(P2)*Y(P1)
       Area     = Area + DArea

    End Do

   !Part 9:
    IF((NEdgeOfCell(Cell)==3 .AND. hasDeadEdge .AND. Area < 0.) .OR. ((NEdgeOfCell(Cell)==4 .OR. (.NOT. hasDeadEdge)) .AND. Area <= 0.))Then
     NegativeVol = 1
     Exit
    EndIF
    
 End Do

!Part 10:
 X(Dead) = TmpX
 Y(Dead) = TmpY
 
!*********************************************************************************************
 End
!###########################################################################################