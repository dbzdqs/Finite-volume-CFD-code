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
 Subroutine Destroy_Tetrahdra(Dim,Heir,IDS,DeadF,LivF,DedF,N1,N2,Cell,NConectCell,IConectCell,&
                             NRemovCell,IRemovCell,NRemovFace,IRemovFace,IFace_Cell,NFace_Cell,ModifiedCells,NModifiedCells,NLivF,NP,NF,NC)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NRemovCell,NRemovFace,DedF,N1,Cell,Cell1,Cell2,Heir,j,I,S,NModifiedCells,NLivF,NF,NC,NP,temp
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:100)::IRemovFace
 Integer,Dimension(1:100)::IRemovCell
 Integer,Dimension(1:8)::DeadF
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:100,1:Dim)::IConectCell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:8)::LivF,N2
 Integer,Dimension(1:100)::ModifiedCells
 Logical::Cellexists,exists
!*********************************************************************************************
!Part 1:
 Cell1 = IDS( N1,DedF )
 Cell2 = IDS( N2(1),LivF(1))

!Part 2:
 Do I=1,NLivF
    If(Cell1==0 .and. N2(I)==1)Then
     temp = IDS(3,LivF(I))
     IDS(3,LivF(I))=IDS(5,LivF(I))
     IDS(5,LivF(I))=temp
     IDS(1,LivF(I))=IDS(2,LivF(I))
     N2(I)=2
    Endif
    IDS(N2(I),LivF(I)) = Cell1
 EndDo
 
!Part 3:
 If(Cell1/=0)Then
  Do J=1,NFace_Cell(Cell1)
     IF( IFace_Cell(J,Cell1) == DedF )Then
      Do I=1,NLivF
         If(I==1)Then
          IFace_Cell(J,Cell1) =  LivF(I) 
         Else
          NFace_Cell(Cell1) = NFace_Cell(Cell1) + 1
          IFace_Cell(NFace_Cell(Cell1),Cell1) =  LivF(I)
         Endif        
      EndDo
     EndIf
  EndDo
    
 !Part 4:
  Cellexists=.FALSE.
  Do S=1,NModifiedCells
     If(ModifiedCells(S)==Cell1)Then
      Cellexists = .TRUE.
      exit
     Endif
  EndDo
  If(.NOT. Cellexists)Then   
   NModifiedCells  =  NModifiedCells + 1
   ModifiedCells(NModifiedCells) = Cell1
  Endif
 Endif
  
!Part 5:
 exists = .FALSE.
 Do J=1,NRemovFace
    If(IRemovFace(J)==DedF)Then
     exists = .TRUE.
     exit
    Endif    
 EndDo
 If(.NOT. exists)Then
  NRemovFace  =  NRemovFace + 1
  IRemovFace(NRemovFace) = DedF
 Endif
 
!Part 6:
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
      NRemovFace = NRemovFace + 1
      IRemovFace(NRemovFace) = DeadF(J)
     Endif
            
    Endif
 EndDo
    
!Part 7:
 exists = .FALSE.
 Do I=1,NRemovCell
    If(IRemovCell(I)==Cell2)Then
     exists = .TRUE.
     exit
    Endif    
 EndDo
 If(.NOT. exists)Then
  NRemovCell  = NRemovCell  +  1
  IRemovCell(NRemovCell) = Cell2
 Endif

!Part 8:
 If(Cell1/=0)Then
  NConectCell(Heir)  =  NConectCell(Heir)  +  1
  IConectCell( NConectCell(Heir) ,Heir) = CELL1
 Endif

!*********************************************************************************************
 End
!###########################################################################################