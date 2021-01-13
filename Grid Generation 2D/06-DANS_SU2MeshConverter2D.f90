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
!// The main objective of this code is to convert SU2 formats to Cell-based or Edge-based data         //!
!// structure. SU2 meshes are used in SU2 open source code.                                            //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************

Program Main_SU2MeshConverter2D
Implicit None
!*********************************************************************************************
!===========================================================================================
!Define the maximum dimensions for the arrays
Integer,Parameter::Dim=99000,DimIDS=4
!===========================================================================================
Integer::MeshDim !Mesh Dimension (2D=2 , 3D=3)
Integer::NC !Number of Cells of mesh
Integer::NP  !Number of Existing Points
Integer::NBoundCrv !Number of Boundary Curves
Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
Integer::NR   !Number Of Regions of mesh
Integer::i
 Integer::Sum                    !Sum of the Number of  Face of each Regions
Character*100,Dimension(1:10):: BCName
Integer,Dimension(1:4,1:Dim)::Corn !Corners point index of Each Cell 
Integer,Dimension(1:4,1:Dim)::Neib !Neighboring Cell Index of Each Elements 
Integer,Dimension(1:4,1:Dim)::NewCorn
Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
Integer,Dimension(1:Dim,1:2)::BFacPt  !Boundary Edge Forming Point
Integer,Dimension(1:100)::NFacCrv !Number of Edges Belong to each Curves
Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh
Character*100,Dimension(1:100)::BCTitle !Boundary Curve Titles
Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
Integer,Dimension(1:100)::BCType !Boundary condition index
Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
Integer,Dimension(1:100)::CellType !Cell Types
!********************************************* main ****************************************
!Part 1: The subroutine for reading the *.SU2 mesh file
 Call Read_2DMeshSU2(Dim,MeshDim,NC,Corn,NP,X,Y,NBoundCrv,BCName,NFacCrv,BFacPt,CellType)
    
!Part 2: Finding the neighbors of each elements
 Call FindNeib2D(Dim,NC,Corn,CellType,Neib)

!Part 3: Write the *.SU2 data in *.plt file 
 Call WritePlaneMesh_cgid_plt(Dim,NP,NC,Corn,X,Y,Z)
 Call WriteBoundCrv_cgid_plt(Dim,NP,NBoundCrv,NFacCrv,BFacPt,X,Y,Z)
    
!Part 4: Converting Cell-Based Mesh to the Edge-Based Mesh and finding the "IDS" from "Corn" and "Neib"
 Call Cell_to_Edge_Hybrid2D(Dim,NC,Corn,Neib,CellType,NF,IDS)
    
!Part 5: Organising the Faces according to their regions
 Call OrganiseIDS(Dim,NBoundCrv,NFacCrv,BFacPt,NF,IDS)
    
!Part 6: Calculating and defining some variables before post-processing for Edge-Based Mesh
!Define Z for 2D
 Z=0
!The number of region is one more than the number of the boundaries (considering interior region)
 NR = NBoundCrv + 1
    
!Define number of faces for each region
 Sum=0
 Do i=1,NR-1
    NFR(i) = NFacCrv(i)
    BCTitle(i) = BCName(i)
    BCType(i) = 2
    Sum = Sum + NFR(i)
 End Do
!Define the interior region
 BCType(i) = 1
 NFR(i) = NF - Sum
 BCTitle(i) = "Interior"
!In 2D all faces are a simple line (edge)
 FaceType(1:NF) = 2
    
!Part 7: Write the Edge-Based Mesh
 Call WriteMesh_gid(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim)
 
End Program Main_SU2MeshConverter2D