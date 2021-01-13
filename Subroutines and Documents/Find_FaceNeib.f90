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
Subroutine Find_FaceNeib(Dim,NF,NC,Dead,Heir,Cell,IDS,NFace_Cell,IFace_Cell,DeadF,LivF,DedF,N1,N2,RemainFaceCount,FaceType,NLivF)
 Implicit None
!*********************************************************************************************
 Intent(In   )::NF,NC,Dead,Heir,Cell,IDS,IFace_Cell,NFace_Cell,FaceType
 Intent(Out  )::DeadF,LivF,DedF,N1,N2,RemainFaceCount,NLivF
 
 Integer::Dim,NF,NC,I,J,K,Face,DedF,Cell,Dead,Heir,P1,P2,P3,P4,NDeadF,NLivF,N1,RemainFaceCount,Point
 Integer,Dimension(1:8)::LivF,N2
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
 Integer,Dimension(1:Dim)::FaceType
 Integer,Dimension(1:8)::DeadF
 LOGICAL::LFace,DFace,exists
!*********************************************************************************************
!Part 1:
 K        = 0
 DeadF(:) = 0
 DedF     = 0
 LivF(:)  = 0
 NLivF    = 0
 NDeadF   = 0

!Part 2:
 RemainFaceCount = NFace_Cell(Cell)   
 
!Part 3:
 Do I=1,NFace_Cell(Cell)
    Face = IFace_Cell(I,Cell)
    P1   = IDS(3,Face)
    P2   = IDS(4,Face)
    P3   = IDS(5,Face)
    If(FaceType(Face)==4)Then
     P4  = IDS(6,Face)
    Else
     P4  = 0
    Endif
    
   !Part 4:
    If( P1/=Dead .and. P2/=Dead .and. P3/=Dead .and. (P4==0 .OR. P4/=Dead) )Then
     NLivF = NLivF + 1
     LivF(NLivF) = Face
    Endif
    
   !Part 5:
    If((P1==Dead .OR. P2==Dead .OR. P3==Dead) .AND. (P1==Heir .OR. P2==Heir .OR. P3==Heir) .AND. P4==0)Then
     NDeadF = NDeadF + 1
	 DeadF(NDeadF)   = Face
     RemainFaceCount = RemainFaceCount - 1
    EndIF
   
   !Part 6:
    If( (P1==Dead .or. P2==Dead .or. P3==Dead .OR. (P4/=0 .AND. P4==Dead)) .and. (P1/=Heir .and. P2/=Heir .and. P3/=Heir .and. (P4==0 .OR. P4/=Heir)) .and. (DedF==0))Then
     DedF = Face
    EndIf
    
 End Do
  
!Part 7:
 If(NLivF==2 .and. DedF==0)Then
  DedF  = LivF(2)
  NLivF = 1   
 Endif
 If(DedF/=0)Then
  If(IDS(1,DedF)==Cell) N1 = 2
  If(IDS(2,DedF)==Cell) N1 = 1
 Endif

!Part 8:
 Do J=1,NLivF
    If(IDS(1,LivF(J))==Cell) N2(J) = 1
    If(IDS(2,LivF(J))==Cell) N2(J) = 2
 EndDo  
!*********************************************************************************************
 End
!###########################################################################################

    