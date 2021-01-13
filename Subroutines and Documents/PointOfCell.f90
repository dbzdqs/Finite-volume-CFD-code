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
!// Date: Feb., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine PointOfCell(Dim,NC,NEdgOfCell,InxEdgOfCell,IDS,Corn)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NEdgOfCell,InxEdgOfCell,IDS
 Intent(out  )::Corn
 
 Integer::Dim,NF,NC,ME,NE,I,J,J1,J2,E,E1,E2,E3,P1_E2,P2_E1,P
 Integer,Dimension(1:4, 1:Dim)::IDS
 Integer,Dimension(1:4, 1:Dim)::InxEdgOfCell,Corn
 Integer,Dimension(1:Dim)::NEdgOfCell
!*********************************************************************************************
 !Part 1:
 do I=1,NC
    Do J=1,NEdgOfCell(I)  
        
       Corn(j,I)=0  
       E = InxEdgOfCell(J,I)
       
      !Part 2:
       IF( E>0 )Then
        P = IDS(3,E)
       Else
        P = IDS(4,-E)
       EndIF

       Corn(J,I) = P

    End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################





