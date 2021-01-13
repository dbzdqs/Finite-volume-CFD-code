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
 Subroutine CheckDiagonalEdges(Dim,Dead,Point,IDS,Facetype,NC,NF,IFace_Cell,NFace_Cell,NConectCellToPoint,IConectCellToPoint,InvalidCell)
 Implicit None
!*********************************************************************************************
 INTENT(IN)::NConectCellToPoint,IConectCellToPoint,NFace_Cell,IFace_Cell,Facetype,IDS,Dead,Point
 INTENT(OUT)::InvalidCell

 Integer::Dim,NConectCellToPoint,A,T,Face2,Cell2,NC,NF,J2,P1C2,P2C2,Dead,Point
 Integer,Dimension(1:100)::IConectCellToPoint
 Logical::InvalidCell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:6,1:Dim)::IDS
!*********************************************************************************************
!Part 1:
 InvalidCell = .FALSE.
 Do A=1,NConectCellToPoint
    Cell2 = IConectCellToPoint(A)
    
   !Part 2:
    Do T=1,NFace_Cell(Cell2)
       Face2 = IFace_Cell(T,Cell2)
       
      !Part 3:
       If(Facetype(Face2)==4)Then
        Do J2=3,6
            
           P1C2  = IDS(J2,Face2)
           If((J2+2)<=6)Then
            P2C2 = IDS((J2+2),Face2)
           Else
            P2C2 = IDS((J2-2),Face2)
           Endif
           
          !Part 4:
           If( (Dead==P1C2 .AND. Point==P2C2) .OR. (Dead==P2C2 .AND. Point==P1C2) )Then
            InvalidCell = .TRUE.
            exit
           Endif
           
        EndDo
        If(InvalidCell) exit
       Endif
       
    EndDo
    If(InvalidCell) exit
 EndDo
!*********************************************************************************************
 End
!###########################################################################################
