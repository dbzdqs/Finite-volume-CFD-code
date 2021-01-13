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
 Subroutine Decrease_CellFaces(Dim,Cell,Dead,Heir,IDS,NRemovFace,IRemovFace,DeadF,RemainFaceCount,&
                            NFace_Cell,IFace_Cell,ModifiedCells,NModifiedCells,DedF,N1,NConectCell,IConectCell,NP,NF,NC)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NRemovCell,NRemovFace,DedF,LivF,N1,N2,Cell,Cell1,Cell2,Heir,j,I,Dead,FTemp,S,K,Point,NP,NF,NC
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:100)::IRemovFace
 Integer,Dimension(1:8)::DeadF
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:100,1:Dim)::IConectCell
 Integer::RemainFaceCount
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:100)::ModifiedCells !Cells that their faces was modified
 Integer::NModifiedCells
 Logical::Cellexists,exists,repeated
!*********************************************************************************************
!Part 1:    
 Do J=1,NFace_Cell(Cell)
    Do I=1,8
       if(Deadf(I)/=0 .AND. IFace_Cell(J,Cell)==Deadf(I))Then
        IFace_Cell(J,Cell) = 0
       Endif
    EndDo
 EndDo

!Part 2:
 Call AscSortArray(IFace_Cell(:,Cell),NFace_Cell(Cell))

!Part 3:
 NFace_Cell(Cell) = RemainFaceCount
 
!Part 4:
 Cellexists = .FALSE.
 Do S=1,NModifiedCells
    If(ModifiedCells(S)==Cell)Then
     Cellexists = .TRUE.
    Endif
 EndDo
 If(.NOT. Cellexists)Then
  NModifiedCells  = NModifiedCells + 1
  ModifiedCells(NModifiedCells) = Cell
 Endif
 
!Part 5:
 Do J=1,8
    If(DeadF(J)/=0)Then     
     exists = .FALSE.
     Do I=1,NRemovFace
        If(IRemovFace(I)==DeadF(J))Then
         exists = .TRUE.
         exit
        Endif    
     EndDo
     If(.NOT. exists)Then
      NRemovFace  =  NRemovFace + 1
      IRemovFace(NRemovFace) = DeadF(J)
     Endif
         
    Endif
 EndDo
 
 if(DedF/=0)Then
     
   !Part 6:
    Cell1 = IDS( N1,DedF )
    If(Cell1/=0)Then
     NConectCell(Heir) = NConectCell(Heir) + 1
     IConectCell( NConectCell(Heir) ,Heir) = CELL1
    
    !Part 7:
     Cellexists = .FALSE.
     Do S=1,NModifiedCells
        If(ModifiedCells(S)==Cell1)Then
         Cellexists = .TRUE.
        Endif
     EndDo
     If(.NOT. Cellexists)Then      
      NModifiedCells  =  NModifiedCells + 1
      ModifiedCells(NModifiedCells) = Cell1
     Endif
    EndIf
    
 Endif
!*********************************************************************************************
 End
!###########################################################################################