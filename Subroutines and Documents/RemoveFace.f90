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
 Subroutine RemoveFace(Dim,Face,NP,NC,NF,IDS,NFace_Cell,IFace_Cell,FaceType,NR,NFR,BeginOfReg,NModifiedCells,ModifiedCells)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Face,NP,NC,BeginOfReg,NR
 Intent(InOut)::NF,IDS,NFace_Cell,IFace_Cell,FaceType,ModifiedCells,NModifiedCells

 Integer::Dim,I,J,J1,NF,NC,NP,Face,E,Fac,Cell,NModifiedCells,NR,Region,RegionLastFace
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:100)::ModifiedCells
 Logical::MEexists,NEexists
 Integer,Dimension(1:100)::NFR,BeginOfReg
!*********************************************************************************************
!Part 1:
 Region = 0
 Do I=1,NR
    If(Face>=BeginOfReg(I) .AND. (Face<= (BeginOfReg(I)+NFR(I)-1)))Then
     Region = I
     exit
    Endif
 EndDo
 RegionLastFace = (BeginOfReg(Region) + NFR(Region) - 1)
 
!Part 2:
 MEexists = .FALSE.
 NEexists = .FALSE.
 Do I=1,NModifiedCells
    If( ModifiedCells(I) == IDS(1,RegionLastFace) )  MEexists = .TRUE.
    If(IDS(2,RegionLastFace)/=0)Then
    If( ModifiedCells(I) == IDS(2,RegionLastFace) ) NEexists = .TRUE.
    EndIf
 EndDo
 If(.NOT. MEexists)Then
  NModifiedCells = NModifiedCells + 1
  ModifiedCells(NModifiedCells) = IDS(1,RegionLastFace)
 Endif
 If( (.NOT. NEexists) .AND. (IDS(2,RegionLastFace)/=0) )Then
  NModifiedCells = NModifiedCells + 1
  ModifiedCells(NModifiedCells) = IDS(2,RegionLastFace)
 Endif
 
!Part 3:
 IDS(:,Face)    = IDS(:,RegionLastFace)
 FaceType(Face) = FaceType(RegionLastFace)

!Part 4:
 Do I=1,NModifiedCells
    Cell=ModifiedCells(I)
    If(Cell==0) Cycle
    Do J1=1,NFace_Cell(Cell)
       Fac = IFace_Cell(J1,Cell)
       If( Fac==RegionLastFace )then
	    IFace_Cell(J1,Cell) = Face
	   EndIf
    EndDo
 EndDo

!Part 5:
 NFR(Region) = NFR(Region)-1
 
!*********************************************************************************************
 End
!###########################################################################################