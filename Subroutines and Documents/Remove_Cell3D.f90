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
 Subroutine Remove_Cell3D(Dim,DelCell,NF,NFace_Cell,IFace_Cell,NC,IDS,&
                          NP,ModifiedCells,NModifiedCells)
 Implicit None
!*********************************************************************************************
 Integer::Dim,I,DelCell,NF,NC,NP,J,K,S,NModifiedCells,Face
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:100)::ModifiedCells !Cells that their faces was modified
 Logical::Cellexists
!*********************************************************************************************
!Part 1:
 Do J=1,NFace_Cell(NC)
    Face = IFace_Cell(J,NC)
    If( IDS(1,Face)==NC ) IDS(1,Face) = DelCell 
    If( IDS(2,Face)==NC ) IDS(2,Face) = DelCell 
 EndDo
 
!Part 2:
 NFace_Cell(DelCell) = NFace_Cell(NC)
 Do I=1,NFace_Cell(NC)
    IFace_Cell(I,DelCell) = IFace_Cell(I,NC)
 End Do
 
!Part 3:
 Cellexists = .FALSE.
 Do S=1,NModifiedCells
    If(ModifiedCells(S)==DelCell)Then
     Cellexists = .TRUE.
     exit
    Endif
 EndDo
 If(.NOT. Cellexists)Then
  NModifiedCells = NModifiedCells + 1
  ModifiedCells(NModifiedCells) = DelCell
 Endif
 
!Part 4:
 NC = NC - 1

!*********************************************************************************************
 End
!###########################################################################################