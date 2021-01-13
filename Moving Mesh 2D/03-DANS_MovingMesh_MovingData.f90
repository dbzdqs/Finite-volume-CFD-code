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
 Program Main_AdaptiveMovingMesh
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
 Real(8)::DT !Time step
 Real(8)::Xoo,Yoo
 Real(8)::Tet,Xp,Yp,X0,Y0
 Real(8)::time       !Simulation Time
 Real(8)::Minf !infinit Flow Mach number
 Real(8)::alfa
 Real(8)::Omega
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::DelX,DelY
 Real(8),Dimension(1:100,1:100)::Vel_X,Vel_Y
 Real(8),Dimension(1:100,1:100)::AnulVelz
 Real(8),Dimension(1:100)::Xo1,Yo1
 Real(8),Dimension(1:100)::MoveTime
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)

!Part 2:
 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Umesh.plt")
 Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Umesh.plt")
 
!Part 2:
 Call Read_MovingData2D(Dim,NR,Nstp,Vel_X,Vel_Y,AnulVelz,Xo1,Yo1)

!Part 3:
 Call BoundPointLabeling(Dim,IDS,NR,NFR,BC,NBP,IBP)
  

 Nstp=3
 DT = 0.1
 
 IStp=0
 Do While(IStp<Nstp)
    IStp=IStp+1
    print*,Nstp,IStp
    
   !Part 4
    Call Cal_Bound_Displac(Dim,Vel_X,Vel_Y,AnulVelz,Xo1,Yo1,IDS,IStp,NR,NFR,BC,DT,X,Y,DelX,DelY)
   
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
    Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Mmesh.plt")
    Call WriteBoundCrv_gid_plt(Dim,NP,NR,NFR,IDS,X,Y)
    
 End Do    

!*********************************************************************************************
 End
!###########################################################################################