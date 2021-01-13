!DDDDDDDDDDDDDDDDDDDDDDDDDDDAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNSSSSSSSSSSSSSSSSSSSSSSSSSS
!//             /////////////       ////////////////    ////////     //////    ////////////////        //!
!//             /////////////      ////////////////    //////////   //////    ////////////////         //!
!//            /////    /////     //////    //////    //////////// //////    /////                     //!
!//           /////    //////    ////////////////    ///////////////////    ////////////////           //!
!//          /////    //////    ////////////////    ////// ////////////               /////            //!
!//         ///////////////    //////    //////    //////   //////////    ////////////////             //!
!//       ///////////////     //////    //////    //////     ////////    ////////////////              //!
!//          Developer            Assistant    in      Numerical             Sciences                  //!
!//----------------------------------------------------------------------------------------------------//!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace Department, Amirkabir University of Technology           //!
!// Chief Developer: N. msnkre, Aerospace eng., Amirkabir University of Technology                     //!
!// Date: October, 14, 2013                                                                            //!
!//                                                                                                    //!
!// The Program is Available Through the Website: www.DANS.ir                                          //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                               //!
!//----------------------------------------------------------------------------------------------------//!
!// Description:                                                                                       //!
!// The main object of this code is mesh refinement. By introducing new point to the domain, new       //!
!// elements will be generated then Delaunay criteria tend to be checked for new elements. As mesh     //!
!// quality decrease by moving mesh, this procedure can be used for improvement mesh quality in moving //!
!// mesh methods. Considerable point is that owing to the high computational cost of edge base data    //!
!// structure, this structure should be changed to a cell-based data structure in the first step and   //!
!// then refinement or coarse mesh generation schemes implemented.                                     //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program MeshRefinement2D
 Implicit None
!*********************************************************************************************
!===============================
 Integer,Parameter::Dim=400000
!===============================
 Integer::NC !Number of Cells of mesh
 Integer::NP  !Number of Existing Points
 Integer::NR   !Number Of Regions of mesh
 Integer::I
 Integer::E1  !Index of Edge between P1 and EP1
 Integer::E2   !Index of Edge between P2 and EP1
 Integer::E3   !Index of Edge between P1 and EP2
 Integer::P1,P2,P3  !Local variables defined for substituting Information of Data Structured 
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NStack !Number of Edges in the Stack List !Number of Edges in the Stack List
 Integer::CandCell  !Candidate cell to be refined  !Candidate cell to be refined
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::Xt,Yt !Coordinate of CENTERS OF CELLS
 Integer,Dimension(1:Dim)::Stack !Index of Edges in the Stack List
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 
!Part 3:
 Do I=1,NC

   !Part 4:
    E1 = InxEdgOfCell(1, I)
    E2 = InxEdgOfCell(2, I)
    E3 = InxEdgOfCell(3, I)

   !Part 5:
    IF( E1>0 )Then
     P1 = IDS(3,E1)
     P2 = IDS(4,E1)
    Else
     P1 = IDS(4,-E1)
     P2 = IDS(3,-E1)
    EndIF

   !Part 6:
    IF( E2>0 )Then
     P3 = IDS(4,E2)
    Else
     P3 = IDS(3,-E2)
    EndIF

   !Part 7:
    Xt(I)=( X(P1)+X(P2)+X(P3))/3
    Yt(I)=( Y(P1)+Y(P2)+Y(P3))/3

 End Do

!Part 8:
 NStack = 0
 
!Part 9:
 Do I=1,NC

   !Part 10:
    NP=NP+1
    X(NP)=Xt(I)
    Y(NP)=Yt(I)

   !Part 11:
    CandCell = I
    Call AddPointToMeshEBased2D(Dim,IDS,NP,NF,NC,CandCell,X,Y,InxEdgOfCell,NStack,Stack)
         
 End Do

!Part 12:
 Call Do_DelaunayEBased2D(Dim,IDS,X,Y,InxEdgOfCell,Stack,NStack)

!Part 13:
 Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"FineM.plt")
!*********************************************************************************************
 Write(*,*) 'Done'
 pause
 End 
!###########################################################################################   





















