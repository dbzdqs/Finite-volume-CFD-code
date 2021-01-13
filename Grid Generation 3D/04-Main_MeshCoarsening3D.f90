!###########################################################################################
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Program Name: MeshCoarsening3D                                                       //!
!// A Program To contraction the suitable Edges in a 3D Mesh                             //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: August, 20, 2017                                                               //!
!//                                                                                      //!
!// N. msnkre , A. Hematizadeh , K. Safari                                               //!
!// Iran,Tehran                                                                          //!
!// OpenMesh@chmail.ir                                                                   //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-CommeFPial Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!###########################################################################################
 Program MeshCoarsening3D
 Implicit None
!*********************************************************************************************
!==================================
 Integer,Parameter::Dim = 200000
!==================================
 Integer::I,J,K,S,A
 Integer::NP
 Integer::NF
 Integer::NC
 Integer::NR
 Integer::Heir
 Integer::Dead
 Integer::NDeformCell
 Integer::Cell
 Integer::temp
 Integer::NContractedEdges
 Integer::counter
 Integer::NConectPoint
 Integer::NConectCell_Loc
 Integer::Point
 Integer::Face,P1,P2
 Real(8)::TmpX,TmpY,TmpZ
 Integer::Beta
 Real(8)::CF
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh
 Integer,Dimension(1:100)::NFR !Number of  Face of each Region
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:100)::BeginOfReg
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:Dim)::NFace_Cell !Number of Faces Forming a Cell
 Integer,Dimension(1:Dim)::NPoint_Cell
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:8,1:Dim)::IPoint_Cell
 Integer,Dimension(1:10,1:Dim)::IFace_Cell !Index of Faces Forming a Cell
 Integer,Dimension(1:100,1:Dim)::IConectCell !Index of Cells connected to each point
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:100)::IConectCell_Loc
 Logical,Dimension(1:Dim)::Bound
 Logical::NegativeVol,InvalidCell,ApexChange
 Integer,Dimension(1:100)::IDeformCell
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
!***************************************** Main ********************************************    
!Part 1:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"Fmesh.plt")
 
!Part 2:
 CF      = 2
 counter = 0
 
!Part 3:
100 Call BoundPointLabeling3D(Dim,Bound,NF,NP,IDS,FaceType)
 Call FaceOfCellV2(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 Call PointOfCell3DV2(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell,NPoint_Cell,IPoint_Cell,FaceType)    
 Call ConectedCell3D(Dim,NC,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NP)
 
 BeginOfReg(:) = 0
 temp=1
 Do I=1,NR
    BeginOfReg(I)=temp
    temp=temp+NFR(I)
 EndDo
 
 Call DetectRemovableEdges3D(Dim,CF,NP,NF,NC,IDS,FaceType,NFR,BeginOfReg,X,Y,Z,Bound,NConectCell,IConectCell,NFace_Cell,IFace_Cell,NContractedEdges,IContractedEdges)
 
!Part 4:
 Do I=1,NContractedEdges
 
   !Part 5:
    Point = IContractedEdges(1,I)
    Dead  = IContractedEdges(2,I)
    Heir  = 0
      
   !Part 6:
    If( NConectCell(Dead)==0 .OR. NConectCell(Dead) > 200)    Cycle
     
   !Part 7:
    NConectCell_Loc = NConectCell(Dead)
    Do J=1,NConectCell_Loc
       IConectCell_Loc(J) = IConectCell(J,Dead)
    End Do
             
   !Part 8:       
    If(Point/=0)Then
     InvalidCell = .FALSE.
     Call CheckDiagonalEdges(Dim,Dead,Point,IDS,Facetype,NC,NF,IFace_Cell,NFace_Cell,NConectCell(Point),IConectCell(:,Point),InvalidCell)
     If(InvalidCell) Cycle
                            
    !Part 9:
     TmpX    = X(Dead)  ; TmpY    = Y(Dead)  ; TmpZ    = Z(Dead)
     X(Dead) = X(Point) ; Y(Dead) = Y(Point) ; Z(Dead) = Z(Point)
     Call VolumeCheck3D(Dim,NF,NC,NP,IDS,X,Y,Z,FaceType,NConectCell_Loc,IConectCell_Loc,NFace_Cell,IFace_Cell,NegativeVol)
     X(Dead) = TmpX ; Y(Dead) = TmpY ; Z(Dead) = TmpZ   
     
     If(.NOT. NegativeVol)Then
                            
     !Part 10:
      Call DeformableCell(Dim,Dead,Point,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NDeformCell,IDeformCell,NP,NC)
               
     !Part 11:
      InvalidCell = .FALSE.
      Call CellValidation(Dim,Dead,Point,NC,NF,InvalidCell,NDeformCell,IDeformCell,NFace_Cell,IFace_Cell,FaceType,IDS)
      If(InvalidCell) Cycle   
      Heir = Point
               
     Endif
    Endif
   
    If(Heir==0) Cycle    
    Print *,Dead
    
   !Part 12:
    Call EdgeContraction3D(Dim,Dead,Heir,NFace_Cell,IFace_Cell,IDS,NF,NC,NR,NFR,BeginOfReg,&
                                NConectCell,IConectCell,NPoint_Cell,IPoint_Cell,NP,FaceType,NDeformCell,IDeformCell)
 
   !Part 13:
    Call ConectedCell3D(Dim,NC,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NP)  
    
    700 continue

 End Do

!Part 14:
 Call RemoveContractedFaces3D(Dim,NF,IDS,FaceType,NR,NFR,BeginOfReg)
 
!Part 15:
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"Cmesh.plt") 
 Print*,'Run is finished','NC=',NC,'NP=',NP,'NF=',NF
 pause
 
 End program MeshCoarsening3D
!###########################################################################################