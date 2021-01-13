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
!// Date: Nov., 15, 2014                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine WriteMesh_gid(Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim)
 Implicit None
!*********************************************************************************************
 Intent(In   )                      ::  Dim,DimIDS,NP,NC,NF,NR,NFR,IDS,X,Y,Z,BCType,BCTitle,CellType,FaceType,MeshDim

 Integer                            ::  Dim,DimIDS,I,J,J1,JJ,NP,NC,NF,NR,SFace,MeshDim,NE,SF
 Integer,Dimension(1:100)           ::  NFR,BCType
 Integer,Dimension(1:DimIDS,1:Dim)  ::  IDS
 Integer,Dimension(1:Dim)           ::  FaceType,CellType 
 Real(8),Dimension(1:Dim)           ::  X,Y,Z
 Character*100,Dimension(1:100)     ::  BCTitle
!*********************************************************************************************
!Part 1:
 Open(1,File='MeshOut.gid')
   
!Part 2:
 Write(1,'(1X,I9,10X,A)') MeshDim  , 'Dimension of Mesh' 

!Part 3:
 Write(1,'(1X,I9,10X,A)') NP       , 'Number of Points'    

!Part 4:
 Write(1,'(1X,I9,10X,A)') NC       , 'Number of Cells'

!Part 5:
 Write(1,'(1X,I9,10X,A)') NF       , 'Number of Faces'

!Part 6:
 Write(1,'(1X,I9,10X,A)') NR       , 'Number of Regions'

!Part 7:
 Write(1,*)'Number of Faces of each Region , BC , Region Name'
 Do J=1,NR
    Write(1,'(15X,I15,20X,I2,3X,A)') NFR(J) , BCType(J) , BCTitle(J)
 End Do 

!Part 8:
 Write(1,*) 'Number of Points of each Face , Left Cell  , Right Cell  ,  Points' 
 Do J=1,NF
    Write(1,'(I15)',Advance='No') FaceType(J)
    Do I=1,FaceType(J)+2            
       Write(1,'(I15)',Advance='No') IDS(I,J)
    End Do
    Write(1,*)
 End Do
	   
!Part 9:
 Write(1,*) 'Coordinates:     X              ,               Y               ,          Z  ' 
 Do J=1,NP
    Write(1,'(10X,ES17.10,5X,ES17.10,5X,ES17.10)') X(J),Y(J),Z(J)
 End Do
	 
!Part 10:
 Write(1,*) '  Cell Type  '
 Do J=1,NC
    Write(1,'(I2)',Advance='No') CellType(J)
    If(Mod(J,20) == 0) Then
     Write(1,*)
    End If
 End Do

 Close(1)
!*********************************************************************************************
 End
!###########################################################################################

