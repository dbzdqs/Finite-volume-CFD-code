!//////////////////////////////////////////////////////////////////////////////////////////!
!// RBF_MovingMesh2D                                                                     //!
!// Date         : Febreury/2/2015                                                       //!
!// Developed by : M. Namvar, Iran, Tehran, OpenFlows@chmail.ir                          //!
!// Version      : V1                                                                    //!
!//                                                                                      //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.Marketxde.ir                       //!
!// It May Be xpied, Modified, And Redistributed For Non-xmmercial Use.                  //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*******************************************************************************************
 Program RBF_MovingMesh3D
 Implicit None
!===============================
 Integer,Parameter::Dim=1000000
!===============================
 Integer::I,J
 real(8)::T1,T2
 
 Integer::NC
 Integer::NP
 Integer::NF
 Integer::NR
 Integer::NStp
 Integer::NBP
 Integer::Move
 Integer::Istp
 Integer::NSBP
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC
 Integer,Dimension(1:100)::BCType
 Integer,Dimension(1:6,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Integer,Dimension(1:Dim)::IBP !Index of Boundary Points
 integer,dimension(1:Dim)::ACTV
 Real(8),Dimension(1:Dim)::X,Y,Z !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8),Dimension(1:Dim)::DelX,DelY,DelZ
 Integer,Dimension(1:Dim)::FaceType !Type of Face (Triangle or Rectangle)
!***************************************** Main ********************************************
!Part 1:
 Call Read_3DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,FaceType,X,Y,Z)

!Part 2:
 Call Write3DMeshSepRgn_gid_plt(Dim,NP,NF,NR,NFR,IDS,X,Y,Z,FaceType,"MeshI.plt")
 Call Write3DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,Z,FaceType)
 Call WriteBoundSurf_gid_plt(Dim,NP,NR,NFR,BC,IDS,X,Y,Z,FaceType)

 !Part 4:
 Istp=0
 Nstp=1
 
 Do While(Istp<Nstp) 
    Istp=Istp+1

   !Part 4:
    Call Read_BDisplacement3D(Dim,NStp,NBP,IBP,DelX,DelY,DelZ)
    
   !Part 5: 
    !!!Call RBF_Moving_Mesh3D(Dim,NBP,NP,IBP,X,Y,Z,DelX,DelY,DelZ)
    
   !For Apply Greedy Algorithm
    Call RBF_GreedyMovingMesh3D(Istp,Dim,NBP,NP,IBP,X,Y,Z,DelX,DelY,DelZ,NSBP,ACTV)
     
   !Part 6:
	Do J=1,NP
       X(J)=X(J)+DelX(J)
       Y(J)=Y(J)+DelY(J)
       Z(J)=Z(J)+DelZ(J)
	End Do

   !Part 7:
    !Call MoveCheckEbased3D(Dim,NF,NC,IDS,X,Y,Z,FaceType,Move)

   !Part 8:
	IF( Move==-1 )Then
	 Do J=1,NP    
        X(J)=X(J)-Delx(J)   
        Y(J)=Y(J)-Dely(J)
        Z(J)=Z(J)-DelZ(J)
     End Do

	 Print*,'Can not Move The mesh'
     pause
	 Exit
	EndIF

   !Part 9:
    Call Write3DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,Z,FaceType)
    Call WriteBoundSurf_gid_plt(Dim,NP,NR,NFR,BC,IDS,X,Y,Z,FaceType)
    
	Print*,'Moving Step:',Istp,'   Predefined Moving Steps:',Nstp

 End Do
!*******************************************************************************************
 pause
 End 
!###########################################################################################
