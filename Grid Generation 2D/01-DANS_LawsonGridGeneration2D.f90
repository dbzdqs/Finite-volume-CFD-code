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
!// Two dimensional Delaunay mesh generation by Lawson method is main object of this code. In this     //!
!// method after introducing new point, location of new point tend to be detected and the triangle that//!
!// this new point is located in will be found. By connecting vertices of triangle to new point, new   //!
!// triangles are generated. Lawson method do not have any warranty for keeping boundaries and         //!
!// preventing the generation of elements outside the desired domain, so to address this problem a     //!
!// blanket and cure-all solution designed .This new method possess considerable advantage over other  //! 
!// methods. Edge base data structure widely used in this code, because it needs lower memory and      //!
!// possess distinguished effectiveness.                                                               //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program Main_LawsonGridGeneration2D
 Implicit None
!===============================
 Integer,Parameter::Dim=80000
!===============================
 Integer::I,J
 Integer::P1,P2,P3 !Local Variables defined for saving Corner Index of each Elements 
 Integer::Inode !!Index of New Node
 Integer::ME !Main Element
 Integer::NE !!Neighboring Element
 Integer::NP  !Number of Existing Points
 Integer::NC !Number of Cells of mesh
 Integer::NBoundCrv !Number of Boundary Curves
 Integer::NStack !Number of Edges in the Stack List !Number of Edges in the Stack List
 Integer::delauny !!Show the Delaunay Criteria 
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NR   !Number Of Regions of mesh
 Integer::MeshDim !Mesh Dimension (2D=2 , 3D=3)
 Integer::DimIDS !first dimension if IDS array (2D=4 , 3D=6)
 Integer::N_Ungrad !Number of Ungrad elements
 Real(8)::Xn,Yn !!Coordinate of the New Introduced point
 Integer,Dimension(1:Dim,1:4)::Corn !Corners point index of Each Cell 
 Integer,Dimension(1:Dim,1:4)::Neib !Neighboring Cell Index of Each Elements 
 integer,dimension(1:Dim,1:2)::Istack !Index of Edges in the Stack List
 integer,dimension(1:Dim,1:2)::BFacPt  !Boundary Edge Forming Point
 Integer,Dimension(1:100)::NFacCrv !Number of Edges Belong to each Curves
 Real(8),Dimension(1:Dim)::X,Y,Z=0.0 !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Integer,Dimension(1:Dim)::Grad !Show the Element is Graded, not Graded or Should be Defined
 Real(8),Dimension(1:Dim)::SF !Size Function
 Integer,Dimension(1:Dim,1:10)::I_Eset !Index of Elements connected to each points
 Integer,Dimension(1:Dim,1:10)::I_Pset  !Number of Points in connected to each points
 Integer,Dimension(1:Dim)::N_Eset   !Number of Elements connected to each points
 Integer,Dimension(1:Dim)::N_Pset   !Number of Points in connected to each points
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BCType !!Boundary condition index
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
 Integer,Dimension(1:Dim)::CellType !Cell Types 
 Character*100,Dimension(1:100)::BCTitle !Boundary Curve Titles
!***************************************** Main ********************************************
!part 1:
 Call Read_2DMeshC(Dim,NP,NC,NBoundCrv,NFacCrv,BFacPt,Corn,Neib,X,Y) 
 
!part 2:
 Call WriteBoundCrv_cgid_plt(Dim,NP,NBoundCrv,NFacCrv,BFacPt,X,Y,Z)
 
!part 3:
 Call Super_Tri(Dim,NP,X,Y,NC,Corn,Neib)

!part 4:
 Do Inode=1,np-3
    print*,inode
   
   !part 5:
    Xn=X(Inode)
	Yn=Y(Inode)

   !part 6:
    call Tri_Contain_Point(Dim,NC,Corn,X,Y,Xn,Yn,ME)

   !part 7:
    call Cons_New_Tri_Lawson(Dim,ME,Inode,NC,Corn,Neib,Nstack,Istack)

   !part 8:
	do while( nstack/=0 )
	
	  !part 9:
	   ME=Istack(nstack,1)
	   NE=Istack(nstack,2)     
	   nstack=nstack-1

      !part 10:
       if( NE/=neib(ME,1) .and. NE/=neib(ME,2) .and. NE/=neib(ME,3) ) cycle
       if( ME/=neib(NE,1) .and. ME/=neib(NE,2) .and. ME/=neib(NE,3) ) cycle

      !part 11:
	   call Del_Check2D(Dim,ME,NE,Corn,Neib,X,Y,delauny)      

      !part 12:
	   if(delauny==-1) call Swapping(Dim,ME,NE,Corn,Neib,Nstack,Istack)

	end do   

 End Do 

!Part 13:
 call BE_Recovery1_2D(Dim,NP,NC,NBoundCrv,NFacCrv,BFacPt,Corn,Neib,X,Y)
pause
!Part 14:
 call Delet_Undesiered_TriPoint(Dim,NC,NP,NBoundCrv,NFacCrv,BFacPt,Corn,Neib,X,Y)

 !Part 15:
  call CellToEdge(Dim,NC,Corn,Neib,NF,IDS)
  call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Omesh.plt")
  print*,"Generating initial mesh has been finished"
 
!part 16:
 call Size_Function(Dim,BFacPt,NFacCrv,NBoundCrv,X,Y,SF)

!Part 17:
 Do i=1,NC
    Grad(i)=0
 End do

!Part 18:
 Do i=1,10000000 
   
   !Part 19:  
    Call Insid_Tri_Grade(Dim,NC,X,Y,SF,Corn,Grad,ME)
    if(ME==0)exit

   !Part 20:
    P1 = Corn(ME,1)
    P2 = Corn(ME,2)
	P3 = Corn(ME,3)

   !Part 21:
    Xn = ( X(P1) + X(P2) + X(P3) ) / 3
	Yn = ( Y(P1) + Y(P2) + Y(P3) ) / 3
 
   !Part 22:
    NP=NP+1
    X(NP)=Xn
	Y(NP)=Yn     
    SF(NP) = ( SF(P1) + SF(P2) + SF(P3) ) / 3
 
   !part 23:
    call Cons_New_Tri_Lawson(Dim,ME,NP,NC,Corn,Neib,Nstack,Istack)
 
   !part 24:
	Grad(NC-1) = 0
	Grad(NC-2) = 0
	Grad(ME  ) = 0

   !part 25:
	do while( nstack/=0 )

      !part 26:
	   ME=Istack(nstack,1)
	   NE=Istack(nstack,2)     
	   nstack=nstack-1

      !part 27:
       if( NE/=neib(ME,1) .and. NE/=neib(ME,2) .and. NE/=neib(ME,3) ) cycle
       if( ME/=neib(NE,1) .and. ME/=neib(NE,2) .and. ME/=neib(NE,3) ) cycle

      !part 28:
	   call Del_Check2D(Dim,ME,NE,Corn,Neib,X,Y,Delauny)      

      !part 29:
	   if(delauny==-1)then
	    call Swapping(Dim,ME,NE,Corn,Neib,Nstack,Istack)
	    Grad(ME) = 0
	    Grad(NE) = 0
	   Endif

	end do   

   !part 30:
    N_Ungrad=0
    Do J=1,NC
	   If(Grad(J)==-1) N_Ungrad=N_Ungrad+1
    End Do
    if(N_Ungrad==1)exit
    
 End Do
 print*,"mesh refinement has been finished"

!Part 31:
 call CellToEdge(Dim,NC,Corn,Neib,NF,IDS)
 call DetectSepetRegnOfMesh2D(Dim,NF,IDS,NR,NFR)
 print*,"converting mesh to edge based data structured has been finished"
 
!Part 32:
 FaceType=2
 MeshDim=2
 DimIDS=4
 Z=0.0
 
!Part 33:
 call WriteMesh_gid(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim)
 call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Omesh.plt")
 print*,"finished"
 
!*********************************************************************************************
 End 
!###########################################################################################
