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
!// The main objective of this code is to extrude of two dimensional mesh in Z direction. The number   //!
!// of layer and the length of extrusion are determined by user.                                       //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program Extrude_2DMesh
 Implicit None
!*********************************************************************************************
!===============================
 Integer,Parameter::Dim=1900000
!===============================
 Integer::I,J
 Integer::ILayer !A counter for Number of layer
 Integer::NF !Number of Faces Constructing Mesh
 Integer::NR   !Number Of Regions of mesh
 Integer::NP  !Number of Existing Points
 Integer::NC !Number of Cells of mesh
 Integer::Nlayer !Number of layers
 Integer::initNP !Initial Number of Points
 Integer::DimIDS !first dimension if IDS array (2D=4 , 3D=6)
 Integer::MeshDim !Mesh Dimension (2D=2 , 3D=3)
 !Real(8)::Thick
 Real(8)::Lenght !Lenght of extrusion determined by user
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::DX,DY,DZ !differentiation in x ,y,z direction 
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:4,1:Dim)::SIDS !!surface IDS(Information of faces Data Structured)
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:4, 1:Dim)::Corn !Corners point index of Each Cell 
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
 Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
 Character*100,Dimension(1:100)::BCTitle !Boundary Curve Titles
!***************************************** Main ********************************************
!Part 1:
 print*,"enter the nlayer:"
 read*,nlayer
 
!Part 2:
 print*,"enter Lenght for extrusion:"
 read*,Lenght
 
!Part 3:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,SIDS,X,Y)
 initNP = NP
 
!Part 4:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,SIDS,X,Y,"InMeh.plt")
 
!Part 5:
 Call EdgeOfCell(Dim,NF,NC,SIDS,NEdgOfCell,InxEdgOfCell)
 
!Part 6:
 Call PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,SIDS,Corn)
 
!Part 7:
 Call Extrude2DMesh(Dim,Nlayer,NF,NR,NFR,BC,NP,NC,SIDS,Corn,NEdgOfCell,FaceType,IDS)
 
!Part 8:
 Do ILayer=1,NLayer	
     
    DX(ILayer) = 0.0
    DY(ILayer) = 0.0
    DZ(ILayer) = Lenght/Nlayer

    Do J=1,initNP
		
       X( J + ILayer*initNP ) = X(J) + DX(ILayer)
       Y( J + ILayer*initNP ) = Y(J) + DY(ILayer)
       Z( J + ILayer*initNP ) = Z(J) + ILayer*DZ(ILayer)

    End Do

 End do
 
!Part 9:
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"ExMeh.plt")
 
!Part 10:
 DimIDS=6
 MeshDim=3
 Call WriteMesh_gid(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BC,BCTitle,FaceType,FaceType,MeshDim)
!*********************************************************************************************
 Print*,"end of program"
 End
!###########################################################################################











