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
MODULE FaceOfCellV3Header  
IMPLICIT NONE
INTERFACE
SUBROUTINE FaceOfCellV3(NF,NC,IDS,NFace_Cell,IFace_Cell,DirFace_Cell)

 Integer::NF,NC
 Integer,Dimension(:,:)::IDS
 Integer,Dimension(:,:)::IFace_Cell
 Integer,Dimension(:)::NFace_Cell
 Integer,Dimension(:,:)::DirFace_Cell
 
END SUBROUTINE FaceOfCellV3
END INTERFACE   
END MODULE FaceOfCellV3Header
!*********************************************************************************************
 Subroutine FaceOfCellV3(NF,NC,IDS,NFace_Cell,IFace_Cell,DirFace_Cell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::NF,NC,IDS
 Intent(Out  )::NFace_Cell,IFace_Cell
 
 Integer::NC,NF,ME,NE,I
 Integer,Dimension(:,:)::IDS
 Integer,Dimension(:)::NFace_Cell
 Integer,Dimension(:,:)::IFace_Cell
 Integer,Dimension(:,:)::DirFace_Cell
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    NFace_Cell(I)=0
 End Do

!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1, I)
    NE = IDS(2, I)
	
   !Part 4:
    NFace_Cell(ME) = NFace_Cell(ME) + 1
    IFace_Cell(NFace_Cell(ME),ME)=I
    DirFace_Cell(NFace_Cell(ME),ME)=1
 
    IF(NE/=0)Then
	 NFace_Cell(NE) = NFace_Cell(NE) + 1
     IFace_Cell(NFace_Cell(NE),NE)=I
     DirFace_Cell(NFace_Cell(NE),NE)=-1
    EndIF

 End Do
!*********************************************************************************************
 End
!###########################################################################################



