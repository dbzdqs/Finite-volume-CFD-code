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
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine FaceOfCellV2(Dim,NF,NC,IDS,NFace_Cell,IFace_Cell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::NF,NC,IDS
 Intent(Out  )::NFace_Cell,IFace_Cell
 
 Integer::Dim,NF,NC,ME,NE,I,J,E1,E2,E3,P1_E2,P2_E1,J1,J2,E,P,N
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:10,1:Dim)::IFace_Cell
 Integer,Dimension(1:Dim)::NFace_Cell
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    NFace_Cell(I) = 0
 End Do

!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1, I)
    NE = IDS(2, I)
	
   !Part 4:
    IF (ME/=0) THEN
	 NFace_Cell(ME) = NFace_Cell(ME) + 1
     IFace_Cell(NFace_Cell(ME),ME)  =  I
    END IF
 
   !Part 5:
    IF(NE/=0)Then
	 NFace_Cell(NE) = NFace_Cell(NE) + 1
     IFace_Cell(NFace_Cell(NE),NE) =  I
    EndIF

 End Do
!*********************************************************************************************
 End
!###########################################################################################



