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
!// The main object of this code is checking Delaunay criteria in mesh elements. As Delaunay criteria  //!
!// maximize the minimum angle of mesh elements, it can be considered as a method of improvement mesh  //!
!// quality. In addition to moving mesh methods and mesh refinement procedure, it can be used for      //!
!// upgrading mesh quality. Considerable point is that owing to the high computational cost of edge    //!
!// base data structure, this structure should be changed to a cell-based data structure in the first  //!
!// step and then refinement or coarse mesh generation schemes implemented.                            //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program MeshImprovmentDelaunay2D
 Implicit None
!*********************************************************************************************
!===============================
 Integer,Parameter::Dim=200000
!===============================
Integer::I
Integer::E1   !Index of Edge between P1 and EP1
Integer::E2   !Index of Edge between P2 and EP1
Integer::E3   !Index of Edge between P1 and EP2
Integer::E4   !Index of Edge between P2 and EP2
Integer::EI   !Selected Edge’s Index In IDS
Integer::P1,P2,P3,P4  !Local variables defined for substituting Information of Data Structured 

Integer::NC !Number of Cells of mesh
Integer::NP  !Number of Existing Points
Integer::NR   !Number Of Regions of mesh
Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
Integer::NStack !Number of Edges in the Stack List !Number of Edges in the Stack List
Integer::CandCell  !Candidate cell to be refined  !Candidate cell to be refined
Integer::DeL !Is Delauny or not
Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
Integer,Dimension(1:Dim)::Stack !Index of Edges in the Stack List
Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
Integer,Dimension(1:100)::BC  !Boundary Condition index
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)

!Part 3:
 NStack = 0
 Do I=1,NF
    NStack = NStack + 1
    Stack(NStack) = I
 End do

!Part 4:
 Do While( NStack/=0 )

   !Part 5:
    EI = Abs( Stack(NStack) )
    NStack = NStack - 1
    IF(IDS(2,EI)==0)Cycle

   !Part 6:
    P1 = IDS(3, EI)
    P2 = IDS(4, EI)

   !Part 7:
    Call Find_PointEdge(Dim,EI,IDS,InxEdgOfCell,P3,P4,E1,E2,E3,E4)

   !Part 8:
    Call DelCeckEBased2D(Dim,P1,P2,P3,P4,X,Y,DeL)

   !Part 9:
    IF(DeL==-1) Call SwapEBased2D(Dim,EI,E1,E2,E3,E4,P1,P2,P3,P4,IDS,InxEdgOfCell,Stack,NStack)

 End do

!Part 10:
 Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"DlnyM.plt")
!*********************************************************************************************
 Write(*,*) 'Done'
 pause
 End 
!###########################################################################################   











