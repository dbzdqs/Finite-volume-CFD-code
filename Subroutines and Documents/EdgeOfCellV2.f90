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
!// Date: Aug., 30, 2015                                                                   //!
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
 Subroutine EdgeOfCellV2(Dim,NF,NC,IDS,NEdgeOfCell,InxEdgeOfCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS
 Intent(Out  )::NEdgeOfCell,InxEdgeOfCell

 Integer::Dim,I,NF,NC,ME,NE
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::InxEdgeOfCell
 Integer,Dimension(1:Dim)::NEdgeOfCell
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    NEdgeOfCell(I)   = 0
    InxEdgeOfCell(:,I) = 0
 End Do

!Part 2: 
 Do I=1,NF
     
   !Part 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
 
   !Part 4:
    NEdgeOfCell(ME) = NEdgeOfCell(ME) + 1
    InxEdgeOfCell(NEdgeOfCell(ME) , ME) = I
 
   !Part 5:
    IF(NE/=0)Then
     NEdgeOfCell(NE) = NEdgeOfCell(NE) + 1
     InxEdgeOfCell(NEdgeOfCell(NE) , NE) = I
    EndIF
    
 End Do
!*********************************************************************************************
 End
!###########################################################################################
