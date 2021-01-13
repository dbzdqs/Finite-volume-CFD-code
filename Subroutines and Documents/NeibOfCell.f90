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
 Subroutine NeibOfCell(Dim,NC,NEdgeOfCell,InxEdgeOfCell,IDS,Neib)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,InxEdgeOfCell,NEdgeOfCell,IDS
 Intent(out  )::Neib
 
 Integer::Dim,NC,I,J,E,N
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::InxEdgeOfCell,Neib
 Integer,Dimension(1:Dim)::NEdgeOfCell
!*********************************************************************************************
!Part 1:
 Do I=1,NC
    Do J=1,NEdgeOfCell(I)
      
      !Part 2:
       Neib(J,I) = 0
 
      !Part 3:
       E = InxEdgeOfCell(J,I)

      !Part 4:
       IF( E>0 )Then
        N = IDS(1,E)
       Else
        N = IDS(2,-E)
       EndIF

      !Part 5:
       Neib(J,I) = N

    End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################





