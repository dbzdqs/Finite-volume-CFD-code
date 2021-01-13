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
!// Date: Feb., 10, 2018                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine FindBoundApexVertices(Dim,NF,NP,IDS,FaceType,NApexVertices,IApexVertices,Bound,X,Y,Z)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NF,NP,I,J,NApexVertices,Face,Point,NBoundaryFaces
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType,IApexVertices,NConnectedFaces
 Logical,Dimension(1:Dim)::Bound
 Integer,Dimension(1:100,1:Dim)::IConnectedFaces
 Integer,Dimension(1:100)::IBoundaryFaces
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:100)::Ux,Uy,Uz,Vx,Vy,Vz
 Real(8)::MaxAngleDiff
!*********************************************************************************************
!Part 1:
 Call FindConnectedFaces(Dim,NF,NP,IDS,FaceType,IConnectedFaces,NConnectedFaces)
 NApexVertices = 0
 Do I=1,NP  
    Point = I
    
   !Part 2:
    If(Bound(Point))Then
     NBoundaryFaces = 0
     Do J=1,NConnectedFaces(Point) 
        Face = IConnectedFaces(J,Point)
        If( IDS(2,Face)==0)Then
         NBoundaryFaces = NBoundaryFaces  +  1
         IBoundaryFaces(NBoundaryFaces) = Face
        EndIF
     End Do
       
    !Part 3:
     Call CalcCrossVectors(Dim,X,Y,Z,Point,FaceType,NP,NF,IDS,NBoundaryFaces,IBoundaryFaces,Ux,Uy,Uz,Vx,Vy,Vz)        
     Call CalcMaxAngleDiff(Dim,NP,X,Y,Z,Point,NBoundaryFaces,Ux,Uy,Uz,Vx,Vy,Vz,MaxAngleDiff)         
     
    !Part 4:
     If(MaxAngleDiff>30)Then
      NApexVertices = NApexVertices   +  1
      IApexVertices(NApexVertices) = Point
     Endif
          
    Endif 
 EndDo
!*********************************************************************************************
 End
!###########################################################################################