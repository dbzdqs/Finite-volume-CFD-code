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
 Subroutine CalcCrossVectors(Dim,X,Y,Z,Point,FaceType,NP,NF,IDS,NBoundaryFaces,IBoundaryFaces,Ux,Uy,Uz,Vx,Vy,Vz)
 Implicit None
!*********************************************************************************************
 Integer::Dim,NP,NF,NBoundaryFaces,I,J,Face,Point,PN,PNV1,PNV2
 Real(8),Dimension(1:Dim)::X,Y,Z
 Integer,Dimension(1:100)::IBoundaryFaces
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
 Real(8),Dimension(1:100)::Ux,Uy,Uz,Vx,Vy,Vz
 Integer,Dimension(1:4)::FacePoints
!*********************************************************************************************
!Part 1:
 Do I=1,NBoundaryFaces
    Face = IBoundaryFaces(I)
    PN   = 0
     
   !Part 2:
    FacePoints(1:FaceType(Face)) = IDS(3:(2+FaceType(Face)),Face)
    Do J=1,FaceType(Face)
       If(Point==FacePoints(J))Then
        PN = J
        exit
       Endif
    EndDo
    
   !Part 3:
    PNV1 = FacePoints(Mod(PN,FaceType(Face)) + 1)
    PNV2 = FacePoints(Mod((PN+FaceType(Face)-2),FaceType(Face))+1) 
    Do J=1,FaceType(Face)
       If(PNV1==FacePoints(J))Then
        PNV1 = J
        exit
       Endif
    EndDo
    Do J=1,FaceType(Face)
       If(PNV2==FacePoints(J))Then
        PNV2 = J
        exit
       Endif
    EndDo
     
   !Part 4:
    Ux(I) = X(FacePoints(PNV1))
    Uy(I) = Y(FacePoints(PNV1))
    Uz(I) = Z(FacePoints(PNV1))
 
    Vx(I) = X(FacePoints(PNV2))
    Vy(I) = Y(FacePoints(PNV2))
    Vz(I) = Z(FacePoints(PNV2))
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################