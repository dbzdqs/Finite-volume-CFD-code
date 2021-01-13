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
!//                                                                                                    //!
!////////////////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************************
 Program MainMovingMesh_MovingEquation
 Implicit None
!===============================
 Integer,Parameter::Dim=100000
!===============================
 
 Integer::J,I,P,IStp 
 Integer::NC !Number of Cells of mesh 
 Integer::NP  !Number of Existing Points 
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh  
 Integer::NR   !Number Of Regions of mesh 
 Integer::NStp 
 Integer::NBP !Number of Boundary Points 
 Integer::Move 
 Integer::Test_Case 
  Integer::NF1,NF2 !Index of 1st and last Non-Boundary Faces
 Integer::NFW1,NFW2 !Index of 1st and last Faces on Wall Boundary 
 Integer::NFF1,NFF2 !Index of 1st and Last Faces on Far-Field Boundary
 Integer::NFI1,NFI2 !Index of 1st and Last Faces on Inflow Boundary 
 Integer::NFS1,NFS2 !Index of 1st ans Last Faces on Symmetry Boundary
 Integer::NFO1,NFO2 !Index of 1st and Last Faces on Inflow Boundary
 Integer::NFIF1,NFIF2 !Index of 1st and last Faces on InterFsce Boundary
 Real(8)::DT !Time step
 Real(8)::Xoo,Yoo
 Real(8)::Tet,Xp,Yp,X0,Y0
 Real(8)::time       !Simulation Time
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::alfa
 Real(8)::Omega
 Real(8)::Total_Time
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::DelX,DelY
 Real(8),Dimension(1:100,1:100)::Vel_X,Vel_Y,Vel_W
 Real(8),Dimension(1:100)::Xo,Yo
 Real(8),Dimension(1:100)::MoveTime
 Integer,Dimension(1:100,1:2)::PtRegion
 Integer,Dimension(1:Dim)::PtIndx
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Umesh.plt")
 Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Umesh.plt")
 
 Call MeshBC(Dim,NR,NFR,BC,IDS,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,NFI1,NFI2,NFS1,NFS2,NFO1,NFO2,NFIF1,NFIF2)
 
!Part 3:
 Call BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)
  
 Time=0.0
 Total_Time=2.2
 IStp=0
 Do While( Time < Total_Time )
    IStp=IStp+1
    print*,IStp,Time,Total_Time
    
   !Part 4
    dT=.1 !time/Nstp
    Test_Case=8
    Minf=0.2
    alfa=0.0
    Omega=0.0
    Call Defin_BoundPoint_Displac(Dim,NFW1,NFW2,IDS,X,Y,Time,dT,Test_Case,Minf,alfa,Omega,Delx,Dely)
    Time=Time+dT
    
   !Part 5:
    Call RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)

   !Part 6:
	Do J=1,NP
       X(J)=X(J)+DelX(J)
       Y(J)=Y(J)+DelY(J)
    End Do

   !Part 7:
    Call MoveCheckEbased2D(Dim,NF,NC,IDS,X,Y,Move)
    if(Move==-1)then
     print*,'negative'
     pause
    endif
  
   !Part 9:
    Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Umesh.plt")
    
 End Do    

!*********************************************************************************************
 End
!###########################################################################################