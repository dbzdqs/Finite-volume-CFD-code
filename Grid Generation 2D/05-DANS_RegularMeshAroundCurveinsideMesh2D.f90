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
!// The main object of this code is the generation of structure and high-density mesh in specific      //!
!// regions. In this code boundary layer mesh generation algorithm used. In this method firstly        //!
!// specific region considered for generating high-density and structure mesh, will be decomposed and  //!
!// then proper curve defined .After that boundary layer mesh tend to be generated and finally         //!
!// duplicated faces will be removed.                                                                  //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
Program InsertCurveToMesh2D
Implicit None
!==============================
Integer,Parameter::Dim = 100000
!==============================
Integer::NR   !Number Of Regions of mesh
Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
Integer::NC !Number of Cells of mesh
Integer::NP  !Number of Existing Points
Integer::NBL  !Number of Boundry Layer mesh
Integer::NBP !Number of Boundary Points
Integer::NPBL  !Number of Point in Boundry Layer mesh
Integer::NC_BL !!Number of Cells of Boundry Layer mesh
Integer::NF_BL !Number of Faces of Boundry layer mesh
Integer::NR_BL !Number of Regions of Boundry Layer mesh
Integer::NE,ME !Neighboring and Main Elements
Integer::NCurv !Number of Curve
Integer::MeshDim !Mesh Dimension (2D=2 , 3D=3)
Integer::NPtCurvs !Number Point on the all of Curves
Integer::FaceType !Type of Face (Triangle or Rectangle)
Integer::NedgCurvs !Number of Edges forming Curves
Integer::I,J,J1,JJ,P1,P2,Pt

Integer,Dimension(1:100)::BC  !Boundary Condition index
Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
Integer,Dimension(1:100)::NFR !Number of  Face of each Regions

Integer,Dimension(100)::DistributionType !type of function to calculate Expansion Ration of Boundry Layer mesh
Integer,Dimension(1:100)::NPtCurv  !Number Points of each Curve
Integer,Dimension(1:100)::NEdgCurv !Number of Edges of each Curve
Integer,Dimension(1:100)::Region_BL !1:generat Boundry Layer mesh for the region, 0: Do not
!Integer,Dimension(1:100)::BL_RegionNum
!Integer,Dimension(1:Dim)::IPtCurvs
Integer,Dimension(1:Dim)::CuspPt!Cusp Point index
Integer,Dimension(1:100)::NFR_BL !Number of Faces in Regions in Boundry Layer mesh
Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
Integer,Dimension(1:2,1:Dim)::EdgPt   !Points index constructing Edges
Integer,Dimension(1:25,1:Dim)::BLPt  !Boundry Layer Point(refer to doc for more information)
Integer,Dimension(1:4,1:Dim)::IDS_BL !Information of mesh Data Structured in Boundry Layer
Integer,Dimension(1:2,1:Dim)::IConectedEdg  !Index of Connected Edge to each point
Real(8),Dimension(1:Dim)::DelX,DelY !Displacement  in x and y Direction !Displacement
Real(8),Dimension(1:Dim)::AngleBP !Angle Boundry Potint
Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
Real(8),Dimension(1:Dim)::XBL,YBL  !Coordinate of Points of Boundry Layer mesh
Real(8),Dimension(1:Dim)::BLThick !Boundry Layer Thickness
!Real(8),Dimension(1:2,1:Dim)::Edgslop
Real(8),Dimension(1:2,1:Dim)::EdgInvSlop !Inverse Slop of Edges 
Real(8),Dimension(1:100)::OveralThick   !Overal Thickness of Boundry Layer mesh 
Real(8),Dimension(1:100)::Dis1  !First Layer Distanse to boundary
Real(8),Dimension(1:100)::Ratio !Expansion Ration of Boundry Layer mesh
Real(8),Dimension(1:100)::StrmThikRatio !Streamwise Thickness Expansion Ratio 
Real(8),Dimension(1:100)::Xref  !Coordinate of References using to modify Streamwise Thickness Expansion Ratio
Integer::NStack !Number of Edges in the Stack List !Number of Edges in the Stack List
Integer,Dimension(1:Dim)::Stack !Index of Edges in the Stack List
Integer,Dimension(1:4,1:Dim)::InxEdgOfCell !Index of Edge Of Cell
Integer,Dimension(1:Dim)::NEdgOfCell !Number of Edge Of Cells
Integer,Dimension(1:Dim)::FaceTyp !Type of Face (Triangle or Rectangle)
Integer::DimIDS !first dimension if IDS array (2D=4 , 3D=6)
Character*100,Dimension(1:100)::BCTitle !Boundary Curve Titles
Integer::InsertCurvReg !Insert Regular mesh around Curve inside Mesh
!*********************************************************************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Input.plt")

!Part 2:
 Call Read_BLayerData2D(Dim,NR,Region_BL,Dis1,OveralThick,Ratio,StrmThikRatio,Xref,NBL,DistributionType)

!Part 3:
 Do I=1,NR
    If(Region_BL(I)/=0) InsertCurvReg = I
 EndDo

!Part 4:
 Call SortFacesOfRegion2D(Dim,InsertCurvReg,NR,NFR,IDS)

!Part 5:
 Call CloneRegn2D(Dim,InsertCurvReg,NP,NF,NR,NFR,IDS,X,Y)
 BC(InsertCurvReg+1)               = BC(InsertCurvReg)
 Region_BL(InsertCurvReg+1)        = Region_BL(InsertCurvReg)
 Dis1(InsertCurvReg+1)             = Dis1(InsertCurvReg)
 OveralThick(InsertCurvReg+1)      = OveralThick(InsertCurvReg)
 Ratio(InsertCurvReg+1)            = Ratio(InsertCurvReg)
 StrmThikRatio(InsertCurvReg+1)    = StrmThikRatio(InsertCurvReg)
 Xref(InsertCurvReg+1)             = Xref(InsertCurvReg)
 DistributionType(InsertCurvReg+1) = DistributionType(InsertCurvReg)

!Part 6:
 Call BL_PreProc(Dim,NR,NFR,IDS,NBL,Region_BL,X,Y,NCurv,NedgCurvs,NPtCurvs,NPBL,XBL,YBL,NEdgCurv,&
                 NPtCurv,EdgPt,IConectedEdg,BLPt,Dis1,OveralThick,Ratio,StrmThikRatio,Xref,DistributionType)

!Part 7:
 Call BoundPtAngle(Dim,NPtCurvs,IConectedEdg,EdgPt,XBL,YBL,AngleBP)

!Part 8:
 CuspPt(:) = 0
 do i=1,NPtCurvs
    if(AngleBP(i)>50)then
     CuspPt(i) = 1
     print*,'cusp point:',i,AngleBP(i)
    endif
 end do

!Part 9:
 do i=1,NPtCurvs
    if( IConectedEdg(1,i)==0 .or. IConectedEdg(2,i)==0 )then
     CuspPt(i) = 1
     print*,'cusp point:',i
    endif
 end do

!Part 10:
 Call InvSlopOfEdges(Dim,XBL,YBL,NedgCurvs,NPtCurvs,IConectedEdg,EdgPt,EdgInvSlop)

!Part 11:
 Call BLayerOveralThickness(Dim,OveralThick,StrmThikRatio,Xref,XBL,NPtCurv,NCurv,BLThick)
 
!Part 12:
 Call GenBLayerLastPt2D(Dim,NPtCurvs,NBL,BLThick,EdgInvSlop,NPBL,BLPt,XBL,YBL)

!Part 13:
 Call GenBLayerPt2D(Dim,NPtCurv,NCurv,DistributionType,BLThick,NBL,BLPt,NPBL,XBL,YBL)

!Part 14:
 Call ExtrudeLineMesh(Dim,NCurv,NedgCurvs,NBL,NEdgCurv,NPtCurvs,EdgPt,IConectedEdg,BLPt,CuspPt,NC_BL,NF_BL,NR_BL,NFR_BL,IDS_BL)

!Part 15:
 Call RemoveRegn2D(Dim,1,2,NPBL,NF_BL,NR_BL,NFR_BL,IDS_BL,XBL,YBL)
 Call Write2DMesh_gid_plt(Dim,NPBL,NF_BL,IDS_BL,XBL,YBL,"BLPrt.plt")

!Part 16:
 If(NC/=0) Call MoveMeshForBoundarylayer(Dim,NC,NR,NF,NP,NPtCurvs,NBL,NFR,BC,IDS,BLPt,XBL,YBL,X,Y)
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"AMPrt.plt")

!Part 17:
 Call MergeMeshes2D(Dim,NP,NF,IDS,NR,NFR,NC,NPBL,NF_BL,IDS_BL,NR_BL,NFR_BL,NC_BL,X,Y,XBL,YBL)
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Merge.plt")

!Part 18:
 FaceTyp=2
 MeshDim=2
 DimIDS=4
 Z=0.0
 call WriteMesh_gid(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BC,BCTitle,FaceTyp,FaceTyp,MeshDim)
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Final.plt")
!*********************************************************************************************
 End
!###########################################################################################