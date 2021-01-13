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
 Subroutine CellValidation(Dim,Dead,Point,NC,NF,InvalidCell,NDeformCell,IDeformCell,NFace_Cell,IFace_Cell,FaceType,IDS)
 Implicit None
!*********************************************************************************************
 INTENT(IN)::NDeformCell,IDeformCell,NFace_Cell,IFace_Cell,NC,NF,FaceType,Dead,Point,IDS
 INTENT(OUT)::InvalidCell

 Integer::Dim,J2,K2,S2,T,Face2,NDeformCell,Cell2,RemainFaceCount,NRemovableFace,NC,NF,P1C2,P2C2,P3C2,P4C2,PF3,Dead,Point,Face3
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:100)::IDeformCell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::FaceType
 LOGICAL::InvalidCell,IsRemovable,HasDiffrentPoint
 Integer,Dimension(1:10)::IRemovableFace
!*********************************************************************************************
!Part 1:
 InvalidCell = .FALSE.
 Do J2=1,NDeformCell
     
   !Part 2:  
    Cell2           = IDeformCell(J2)
    RemainFaceCount = NFace_Cell(Cell2)
    NRemovableFace  = 0
    
   !Part 3:  
    Do K2=1,NFace_Cell(Cell2)
        
      !Part 4:  
       Face2 = IFace_Cell(K2,Cell2)
       If(FaceType(Face2)==3)Then
        P1C2 = IDS(3,Face2)
        P2C2 = IDS(4,Face2)
        P3C2 = IDS(5,Face2)
        IF( (P1C2==Dead .or. P2C2==Dead .or. P3C2==Dead) .and. (P1C2==Point .or. P2C2==Point .or. P3C2==Point))Then
         RemainFaceCount = RemainFaceCount  - 1
         NRemovableFace  = NRemovableFace   + 1
         IRemovableFace(NRemovableFace) = Face2
        EndIf
       EndIf
                
    End Do
    
   !Part 5:
    If(RemainFaceCount==3)Then
     InvalidCell = .TRUE.
     return
    EndIf
   
   !Part 6:
    If(RemainFaceCount>2)Then
     InvalidCell = .FALSE.
     
    !Part 7:
     Do K2=1,NFace_Cell(Cell2)
        Face2 = IFace_Cell(K2,Cell2)
         
       !Part 8:
        IsRemovable = .FALSE.
        Do T=1,NRemovableFace
           If(Face2==IRemovableFace(T))Then
            IsRemovable = .TRUE.
            exit
           Endif
        EndDo
        If(IsRemovable) Cycle
        
       !Part 9:
        P1C2 = IDS(3,Face2)
        P2C2 = IDS(4,Face2)
        P3C2 = IDS(5,Face2)
        If(FaceType(Face2)==4)Then
         P4C2 = IDS(6,Face2)
        Else
         P4C2 = 0
        Endif     
        If(P1C2==Dead) P1C2 = Point
        If(P2C2==Dead) P2C2 = Point
        If(P3C2==Dead) P3C2 = Point
        If(P4C2==Dead) P4C2 = Point
        
       !Part 10:
        Do S2=1,NFace_Cell(Cell2)
           Face3 = IFace_Cell(S2,Cell2)
           If(Face3==Face2) Cycle
            
          !Part 11:
           IsRemovable = .FALSE.
           Do T=1,NRemovableFace
              If(Face3==IRemovableFace(T))Then
               IsRemovable = .TRUE.
               exit
              Endif
           EndDo
           If(IsRemovable) Cycle
            
          !Part 12:
           HasDiffrentPoint = .FALSE.
           Do T=3,(2+Facetype(Face3))
              If(IDS(T,Face3)==Dead)Then
               PF3 = Point
              Else
               PF3 = IDS(T,Face3)
              Endif
                            
              If(PF3/=P1C2 .AND. PF3/=P2C2 .AND. PF3/=P3C2 .AND. PF3/=P4C2 )Then
               HasDiffrentPoint = .TRUE.
               exit
              Endif
           EndDo
            
          !Part 13:
           If(.NOT. HasDiffrentPoint)Then
            InvalidCell = .TRUE.
            return
           Endif     
                        
        EndDo
        
       !Part 14:
        If(InvalidCell) return
        
     EndDo
                
     If(InvalidCell) return
                
    Endif    
 EndDo
!*********************************************************************************************
 END 
!###########################################################################################
