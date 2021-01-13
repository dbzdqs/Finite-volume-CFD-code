!###########################################################################################
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Program Name: MeshCoarsening2D                                                       //!
!// A Program to Remove suitable elements of an array of Edges repeatedly in a 2D mesh   //!
!//                                                                                      //!
!// Version: V1                                                                          //!
!// Date: July, 26, 2017                                                                 //!
!//                                                                                      //!
!// N. msnkre , A. Hematizadeh , K. Safari                                               //!
!// Iran,Tehran                                                                          //!
!// OpenMesh@chmail.ir                                                                   //!
!//                                                                                      //!
!// This Program is Available Through the Website: www.MarketCode.ir                     //!
!// It May be Copied, Modified, and Redistributed For Non-Commerical Use.                //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!###########################################################################################
 program MeshCoarsening_2D
 Implicit None
!*********************************************************************************************
!===============================
 Integer,Parameter::Dim=30000  
!===============================
 Integer::I,J,JJ,NP,NF,NFT,NC,EI,NR,NConectCell,NContractedEdge=0,Edge,Heir,Dead,NegetiveVol,P1,P2,E,Edgeexist,NContractedEdges,counter,temp,step
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:100)::BeginOfReg
 Integer,Dimension(1:4,1:Dim)::IEdgeOfCell
 Integer,Dimension(1:Dim)::NEdgOfCell
 Integer,Dimension(1:Dim)::NConectEdge,IDeadPoints
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Integer,Dimension(1:1000)::IConectCell
 Integer,Dimension(1:Dim,1:2)::IContractedEdge
 Integer::NBP !Number of Boundary Points
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Logical::IsInvalidEdge
 Real(8)::CF,Area,DArea,AngleDeg
 Integer,Dimension(1:2,1:Dim)::IContractedEdges
 Integer,Dimension(1:100)::IBoundaryFaces
!***************************************** Main ********************************************    
!Part 1:
 CF = 2.0
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Fmesh.plt")    
 
 step = 1
 BeginOfReg(:) = 0
 temp = 1
 Do I=1,NR
    BeginOfReg(I) = temp
    temp = temp + NFR(I)
 EndDo
 NFT = NF
 
!Part 2:
 Call EdgeOfCellV3(Dim,NF,NC,BeginOfReg,NR,NFR,IDS,NEdgOfCell,IEdgeOfCell)
 
!Part 3:
 Call DetectRemovableEdges(CF,NP,NF,NC,IDS,X,Y,NEdgOfCell,IEdgeOfCell,NContractedEdges,IContractedEdges)
 
!Part 5:
 Call EliminateEdge(Dim,NContractedEdges,IContractedEdges,NEdgOfCell,IEdgeOfCell,BeginOfReg,NR,NFR,NP,NF,NC,IDS,X,Y,NBP,IBP)
 Call RemoveContractedFaces(NF,IDS,NR,NFR,BeginOfReg)
 
!Part 6:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Cmesh.plt")  
 
!Part 7:
 Print*,'Run is finished',NC,NP,NF
 pause
 
 End program MeshCoarsening_2D
!###########################################################################################

