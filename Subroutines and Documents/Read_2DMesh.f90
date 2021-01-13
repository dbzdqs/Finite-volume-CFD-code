!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//       /////////////       ////////////////    ////////     //////    ////////////////  //!
!//       /////////////      ////////////////    //////////   //////    ////////////////   //!
!//      /////    /////     //////    //////    //////////// //////    /////               //!
!//     /////    //////    ////////////////    ///////////////////    ////////////////     //!
!//    /////    //////    ////////////////    ////// ////////////               /////      //!
!//   ///////////////    //////    //////    //////   //////////    ////////////////       //!
!// ///////////////     //////    //////    //////     ////////    ////////////////        //!
!//    Developer            Assistant    in      Numerical             Sciences            //!
!//----------------------------------------------------------------------------------------//!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Supervisor: Dr. h. hdhrnuidn, Aerospace eng. Amirkabir University of Technology      //!
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Read_2DMesh(Dim,NP,NC,NF,NR,NFR,BC,IDS,X,Y)
 Implicit None
!*********************************************************************************************
 Integer::I,J,J1,JJ,SFace,FaceType,MeshDim

 Integer,Intent(In) ::Dim
 Integer,Intent(Out)::NP !number of points
 Integer,Intent(Out)::NC
 Integer,Intent(Out)::NF !Number of Faces Constructing Mesh
 Integer,Intent(Out)::NR
 Integer,Dimension(1:100)    ,Intent(Out)::NFR
 Integer,Dimension(1:100)    ,Intent(Out)::BC
 Integer,Dimension(1:4,1:Dim),Intent(Out)::IDS
 Real(8),Dimension(1:Dim)    ,Intent(Out)::X,Y
!*********************************************************************************************
!Part 1:
 Open(1,File='Mesh.gid')

!Part 2:
 Read(1,*) MeshDim
 IF(MeshDim/=2)Print*,'Please Check the Mesh File. It is not a 2D Mesh'

!Part 3:
 Read(1,*) NP    

!Part 4:
 Read(1,*) NC

!Part 5:
 Read(1,*) NF

!Part 6:
 Read(1,*) NR

!Part 7:
 Read(1,*)
 Do J=1,NR
    Read(1,*) NFR(J) , BC(J)
 End Do 

!Part 8:
 Read(1,*)
 Do J=1,NF  
    Read(1,*) FaceType,IDS(1,J),IDS(2,J),IDS(3,J),IDS(4,J)
 End Do
 
!Part 9:
 Read(1,*)        
 Do J=1,NP
    Read(1,*) X(J),Y(J)
 End Do

 Close(1)
!*********************************************************************************************
 End
!###########################################################################################
