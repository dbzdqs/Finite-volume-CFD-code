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
 Program MainMovingMesh_SomeRegion
 Implicit None
!===============================
 Integer,Parameter::Dim=40000  
!===============================

 Integer::I,k,J,IStp,Pt
 Integer::NC !Number of Cells of mesh
 Integer::NP  !Number of Existing Points
 Integer::NF !Number of Faces Constructing Mesh  !Number of Faces Constructing Mesh 
 Integer::NR   !Number Of Regions of mesh
 Integer::NStp
 Integer::NBP !Number of Boundary Points
 Integer,Dimension(1:100)::NFR !Number of  Face of each Regions
 Integer,Dimension(1:100)::BC  !Boundary Condition index
 Integer,Dimension(1:100)::BCType !Boundary condition index
 Integer,Dimension(1:100)::MoveRegion
 Integer,Dimension(1:4,1:Dim)::IDS !Information of Data Structured (1,i):Left Cell,(2,i):Right Cell,(3:FaceType,i):Forming Point of Face
 Real(8),Dimension(1:Dim)::X,Y !Coordinate of Points Constructing Mesh !Coordinate of Points Constructing Mesh
 Real(8)::Xo,Yo,Tet
 
 Integer,Dimension(1:100,1:2)::PtRegion
 Integer,Dimension(1:Dim)::PtIndx
!***************************************** Main ********************************************
!Part 1:
 Call Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
    

 Call Write2DMeshSepRgn_gid_plt(Dim,NP,NR,NFR,IDS,X,Y,"Umesh.plt")
 Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Umesh.plt")
 
 
 Call RegionsPointData(Dim,NR,NFR,IDS,PtRegion,PtIndx)
 
 
 MoveRegion(1)=0
 MoveRegion(2)=1
 MoveRegion(3)=0
 MoveRegion(4)=0
 MoveRegion(5)=0
 MoveRegion(6)=0

 Xo=0.0 ; Yo=0.0

 Nstp=3
 IStp=0
 Do While(IStp<Nstp)
    IStp=IStp+1
    print*,Nstp,IStp
    
    Do I=1,NR
       IF(MoveRegion(I)==1)Then
        DO J=PtRegion(I,1),PtRegion(I,2)
           Pt=PtIndx(J)
           
           Tet=5.0
           Call RotatePt(Xo,Yo,Tet,X(Pt),Y(Pt))
        END DO
       EndIF
    End Do

    
   !Part 9:
    Call Write2DMesh_gid_plt(Dim,NP,NF,IDS,X,Y,"Umesh.plt")
 
 End Do  
 
!*********************************************************************************************
 End 
!###########################################################################################
