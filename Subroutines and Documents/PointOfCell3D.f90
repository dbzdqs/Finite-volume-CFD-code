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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PointOfCell3D(Dim,FaceType,NFace_Cell,IFace_Cell,IDS,Cell,Corn)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,FaceType,NFace_Cell,IDS,Cell
 Intent(Out  )::Corn
 Intent(InOut)::IFace_Cell
 
 Integer::Dim,I,J,k,P1,P2,NFace,Face,Face1,Face2,Cell,Temp
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:2)::IFace
 Integer,Dimension(1:8)::P
 Integer,Dimension(1:Dim)::FaceType,NFace_Cell
 Integer,Dimension(1:8,1:Dim)::Corn
 Integer,Dimension(1:6,1:Dim)::IFace_Cell
!*********************************************************************************************
!Part 1:
 IF( NFace_Cell(Cell)>4 )THEN
  DO J=1, NFace_Cell(Cell)
     
     !Part 2:
      IF (Facetype(J)==4) THEN
      Temp=IFace_Cell(J,Cell)
      IFace_Cell(J,Cell)=IFace_Cell(1,Cell)
      IFace_Cell(1,Cell)= Temp 
     END IF 
     
  END DO
 END IF
  
!Part 3: 
 P(:) = 0

!Part 4:
 Face = IFace_Cell(1,cell)
 Do I=1,FaceType(Face)
     P(I) = IDS(2+I,Face)
 End Do

 IF( P(4)==0 ) P(4) = P(3)

!Part 5:
 Do k=1,4

    NFace = 0
   
   !Part 6:
    Do I=2,NFace_Cell(Cell)

       Face = IFace_Cell(I,Cell)
 
      !Part 7:
       Do J=1,FaceType(Face)
         
          IF(P(k)==IDS(2+J,Face) )Then
 	       NFace=NFace+1
	       IFace(NFace) = Face
	      EndIF

       End Do

    End Do

!Part 8:
    Face1 = IFace(1)
    Face2 = IFace(2)
 
 !Part 9:
    Do I=1,FaceType(Face1)
       P1 = IDS(2+I,Face1)

       Do J=1,FaceType(Face2)
          P2 = IDS(2+J,Face2)
!Part 10:
          IF(P1==P2 .and. P1/=P(k) ) P(4+k) = P1
       End Do

    End Do

 End Do
 
 !Part 11:
 Corn(:,Cell) = P(:)
   
!*********************************************************************************************
 End
!###########################################################################################



 