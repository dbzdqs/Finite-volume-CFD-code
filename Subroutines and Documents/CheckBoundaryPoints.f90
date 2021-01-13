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
 Subroutine CheckBoundaryPoints(Dim,IDS,NBP,IBP,Dead,Heir,EI,NEdgeOfCell,IEdgeOfCell,IsInvalidEdge)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NBP,IBP,EI,NEdgeOfCell,IEdgeOfCell
 Intent(InOut)::Dead,Heir,IDS
 Intent(Out  )::IsInvalidEdge
    
 Integer::Dim,Dead,Heir,NBP,tmp1,EI,J,I,I3,I4,Cell,CellBoundaryEdgesCount
 Logical::IsInvalidEdge
 Integer,Dimension(1:Dim)::IBP
 Integer,Dimension(1:4,1:Dim)::IDS
 Logical DeadIsBoundary,HeirIsBoundary
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(4,Dim)::IEdgeOfCell
!*********************************************************************************************
!Part 1:
 IsInvalidEdge  = .FALSE.
 DeadIsBoundary = .FALSE.
 HeirIsBoundary = .FALSE.
 
!Part 2:
 Do J=1,NBP
     
   !Part 3:
    IF( IBP(J)==Dead )Then
     DeadIsBoundary = .TRUE.
     exit
    Endif
    
 End Do

!Part 4:
 If( DeadIsBoundary )Then
     
 !Part 5:
  Do J=1,NBP
     IF( IBP(J)==Heir )Then
      HeirIsBoundary = .TRUE.
      exit
     Endif
  End Do
  
 Endif  
    
!Part 6:
 If( DeadIsBoundary .AND. HeirIsBoundary )Then
     
 !Part 7:
  Do I3=1,2
      
    !Part 8:
     Cell = IDS(I3,EI)
     If(Cell/=0)Then
         
     !Part 9:
      CellBoundaryEdgesCount = 0
      Do I4=1,NEdgeOfCell(Cell)
          
        !Part 10:
         if(IDS(2,IEdgeOfCell(I4,Cell))==0)  CellBoundaryEdgesCount = CellBoundaryEdgesCount+1
         
      EndDo
     
     !Part 11:
      If( CellBoundaryEdgesCount > 1 )Then
       IsInvalidEdge = .TRUE.
       return
      Endif
      
     Endif
     
  EndDo
  
!Part 12:
 ElseIf(DeadIsBoundary .AND. (.NOT.(HeirIsBoundary)))Then
  tmp1 = Dead
  Dead = Heir
  Heir = tmp1      
 Endif
 
!*********************************************************************************************
 End
!###########################################################################################