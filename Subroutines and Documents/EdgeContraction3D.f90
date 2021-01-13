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
 Subroutine EdgeContraction3D(Dim,Dead,Heir,NFace_Cell,IFace_Cell,IDS,NF,NC,NR,NFR,BeginOfReg,NConectCell,&
                    IConectCell,NPoint_Cell,IPoint_Cell,NP,FaceType,NDeformCell,IDeformCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim
 Intent(InOut)::Dead,Heir,NFace_Cell,IFace_Cell,IDS,NF,NC,NR,NFR,NConectCell,&
                IConectCell,NPoint_Cell,IPoint_Cell,NP,FaceType,NDeformCell,IDeformCell

 Integer::Dim,I,J,N,k,S,NP,NF,NC,Heir,Dead,N1,F,NDeformCell,Point,NLivF,lastRepeatedLoc,RepeatedCount,NR
 Integer::NRemovCell,NRemovFace,DedF,Cell,Face,DelCell,FTemp,RemainFaceCount,NModifiedCells
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:Dim)::NPoint_Cell
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:8,1:Dim)::IPoint_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:100,1:Dim)::IConectCell
 Integer,Dimension(1:100)::IDeformCell
 Integer,Dimension(1:8)::LivF,N2
 Integer,Dimension(1:100)::IRemovFace
 Integer,Dimension(1:100)::IRemovCell
 Integer,Dimension(1:8)::DeadF
 Logical::repeated,MEexists, NEexists
 Integer,Dimension(1:100)::ModifiedCells
 Integer,Dimension(1:100)::NFR,BeginOfReg
!********************************************************************************************* 
!Part 1:
 NModifiedCells = 0
 NRemovFace     = 0 
 NRemovCell     = 0 
 
!Part 2:
 Do I=1,NDeformCell
    Cell = IDeformCell(I)
    Call Find_FaceNeib(Dim,NF,NC,Dead,Heir,Cell,IDS,NFace_Cell,IFace_Cell,&
            DeadF,LivF,DedF,N1,N2,RemainFaceCount,FaceType,NLivF)
    If(RemainFaceCount<=2)Then
     Call Destroy_Tetrahdra(Dim,Heir,IDS,DeadF,LivF,DedF,N1,N2,Cell,NConectCell,&
                            IConectCell,NRemovCell,IRemovCell,NRemovFace,IRemovFace,IFace_Cell,&
                            NFace_Cell,ModifiedCells,NModifiedCells,NLivF,NP,NF,NC)
    Else
     Call Decrease_CellFaces(Dim,Cell,Dead,Heir,IDS,NRemovFace,IRemovFace,&
                            DeadF,RemainFaceCount,NFace_Cell,IFace_Cell,ModifiedCells,&
                            NModifiedCells,DedF,N1,NConectCell,IConectCell,NP,NF,NC)
    Endif
 EndDo

!Part 3:
 Do I=1,NConectCell(Dead)
    Cell = IConectCell(I,Dead) 
    Do J=1,NFace_Cell(Cell)
        
      !Part 4:
	   Face     = IFace_Cell(J,Cell)
       MEexists = .FALSE.
       NEexists = .FALSE.
       Do S=1,NModifiedCells
          If( ModifiedCells(S) == IDS(1,Face) )  MEexists = .TRUE.
          If(IDS(2,Face)/=0)Then
           If( ModifiedCells(S) == IDS(2,Face) ) NEexists = .TRUE.
          EndIf
       EndDo
       If(.NOT. MEexists)Then
        NModifiedCells = NModifiedCells + 1
        ModifiedCells(NModifiedCells) = IDS(1,Face)
       Endif
       If( (.NOT. NEexists) .AND. (IDS(2,Face)/=0) )Then
        NModifiedCells = NModifiedCells + 1
        ModifiedCells(NModifiedCells) = IDS(2,Face)
       Endif
        
      !Part 5:
       Do k=3,(2+FaceType(Face))
          IF( IDS(k,Face)==Dead ) IDS(k,Face) = Heir
       End Do
        
      !Part 6:
       RepeatedCount   = 0
       lastRepeatedLoc = 0
       Do k=3,(2+FaceType(Face))
          If(IDS(k,Face)==Heir)Then
           RepeatedCount = RepeatedCount + 1
           lastRepeatedLoc = K
          EndIf
       EndDo
       
      !Part 7:
       If(RepeatedCount>1)Then
        Do K=lastRepeatedLoc,(1+FaceType(Face))
           IDS(K,Face) = IDS(K+1,Face)
        EndDo     
        IDS((2+FaceType(Face)),Face) = 0
        FaceType(Face) = FaceType(Face) - 1
       EndIf
        
    EndDo
 EndDo

!Part 8:
 Call AscSortArray(IRemovFace,NRemovFace)
 Do I=1,NRemovFace
    Face = IRemovFace(I)
    Call RemoveFace(Dim,Face,NP,NC,NF,IDS,NFace_Cell,IFace_Cell,FaceType,NR,NFR,BeginOfReg,NModifiedCells,ModifiedCells)
 End Do

!Part 9:
 Call AscSortArray(IRemovCell,NRemovCell)
 Do I=1,NRemovCell
    DelCell = IRemovCell(I)
    Call Remove_Cell3D(Dim,DelCell,NF,NFace_Cell,IFace_Cell,NC,IDS,NP,ModifiedCells,NModifiedCells)
 End Do
 
!Part 10:
 Do I=1,NModifiedCells
    Cell = ModifiedCells(I)
    
   !Part 11:
    If(Cell > NC .OR. Cell==0) Cycle
    
   !Part 12:
    NPoint_Cell(Cell) = 0
    Do J=1,NFace_Cell(Cell)
       FTemp = IFace_Cell(J,Cell)
       Do S=3,2+FaceType(FTemp)
          Point    = IDS(S,FTemp)
          repeated = .FALSE.
          Do K=1,NPoint_Cell(Cell)
             if(IPoint_Cell(K,Cell)==Point)Then
              repeated = .TRUE.
              exit
             Endif
          EndDo
          If(.NOT. repeated)Then
           NPoint_Cell(Cell) = NPoint_Cell(Cell) + 1
           IPoint_Cell(NPoint_Cell(Cell),Cell) = Point
          Endif
       EndDo
    EndDo
    
 EndDo
!*********************************************************************************************
100 End
!###########################################################################################