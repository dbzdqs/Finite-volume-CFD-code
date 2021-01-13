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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine EdgeOfCellV3(Dim,NF,NC,BeginOfReg,NR,NFR,IDS,NEdgeOfCell,IEdgeOfCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS
 Intent(Out  )::NEdgeOfCell,IEdgeOfCell

 Integer::Dim,I,J,NF,NC,ME,NE,NR
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::IEdgeOfCell
 Integer,Dimension(1:Dim)::NEdgeOfCell
 Integer,Dimension(1:100)::NFR,BeginOfReg
!*********************************************************************************************
!Part 1:
 Do I=1,Dim
    NEdgeOfCell(I)   = 0
    IEdgeOfCell(:,I) = 0
 End Do

!Part 2: 
 Do J=1,NR
     Do I=BeginOfReg(J),BeginOfReg(J)+NFR(J)-1
     
       !Part 3:
        ME = IDS(1,I)
        NE = IDS(2,I)
 
       !Part 4:
        NEdgeOfCell(ME) = NEdgeOfCell(ME) + 1
        IEdgeOfCell(NEdgeOfCell(ME) , ME) = I
 
       !Part 5:
        IF(NE/=0)Then
         NEdgeOfCell(NE) = NEdgeOfCell(NE) + 1
         IEdgeOfCell(NEdgeOfCell(NE) , NE) = I
        EndIF
     EndDo
     
 End Do
!*********************************************************************************************
 End
!###########################################################################################
