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
!// The main objective of this code is for extracting Voronoi diagram. By drawing Voronoi diagram in   //! 
!// computational meshes, in spite of decreasing the number of cells, the number of faces will be      //!
!// increased. In this case by reducing the number of cells the computational time decreased and the   //! 
!// accuracy of the solution improved.                                                                 //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program Voronoi2D
 Implicit None
!===========================================================================================
 Integer,Parameter::Dim=50000
!===========================================================================================
 Integer::NP  !Number of Existing Points
 Integer::NC !Number of Cells of mesh
 Integer::NR   !Number Of Regions of mesh
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NewNP !New number of Existing Points for constructing Voronoi Diagram
 Integer::NewNC !New number of Existing Cells for constructing Voronoi Diagram
 Integer::NewNF !New number of Existing faces for constructing Voronoi Diagram
 Integer::MeshDim !Mesh Dimension (2D=2 , 3D=3)
 Integer::NewNR !New number of Regions for constructing Voronoi Diagram
 Integer::i
 Integer::DimIDS !first dimension if IDS array (2D=4 , 3D=6)
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:100)::NewBC !New Boundry conditions for constructing Voronoi Diagram
 Integer,Dimension(1:100)::NewNFR !new Number of  Face of each Regions for constructing Voronoi Diagram
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:4,1:Dim)::NewIDS ! new Information of Data Structured for constructing Voronoi Diagram
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::XCen,YCen !coordinates of Cell Centroid of elements
 Real(8),Dimension(1:Dim)::NewX,NewY !New Coordinate of Points Constructing Voronoi Diagram
 Character*100,Dimension(1:100)::BCTitle !Boundary Curve Titles
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:Dim)::CellType !Cell Types
 !Integer,Dimension(1:Dim)::CellTypeVor
!********************************************* main ****************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
         
!Part 2:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"IMesh.plt")
 
!Part 3:
 Z(1:NP)=0
 MeshDim = 2
 FaceType(1:NF) = 2

!Part 4:
 Call CalCellCentroid(Dim,NF,NC,IDS,X,Y,XCen,YCen)
    
!Part 5:
 Call DualMeshGenerator2D(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y,XCen,YCen,NewNP,NewNC,NewNF,NewNR,NewBC,NewNFR,NewIDS,NewX,NewY)
    
!Part 6:
 DimIDS=4
 BCTitle="..."
 FaceType(1:NewNF) = 2
 CellType(1:NewNC)=-1
 Call WriteMesh_gid(Dim,DimIDS,NewNP,NewNC,NewNF,NewNR,NewNFR,NewIDS,NewX,NewY,Z,NewBC,BCTitle,CellType,FaceType,MeshDim)
     
!Part 7:
 Call Write2DMeshSepRgn_gid_plt(Dim,NewNP,NewNR,NewNFR,NewIDS,NewX,NewY,"OMesh.plt")
!*********************************************************************************************
 End
!########################################################################################### 