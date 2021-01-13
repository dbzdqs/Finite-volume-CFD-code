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
!// Developed by: F. Farhadkhani, Mathmatical, Amirkabir university of Technology          //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine EdgeOfCell(Dim,NF,NC,IDS,NEdgOfCell,InxEdgOfCell)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS
 Intent(Out  )::NEdgOfCell,InxEdgOfCell
 
 Integer::Dim,NF,NC,ME,NE,I,J,E1,E2,E3,P1_E2,P2_E1,J1,J2,E,P,N
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:4,1:Dim)::InxEdgOfCell
 Integer,Dimension(1:Dim)::NEdgOfCell
!*********************************************************************************************

!Part 1:
 Do I=1,NC
    NEdgOfCell(I)=0
 End Do

!Part 2:
 Do I=1,NF

   !Part 3:
    ME = IDS(1,I)
    NE = IDS(2,I)
	
   !Part 4:
	NEdgOfCell(ME) = NEdgOfCell(ME) + 1
    InxEdgOfCell(NEdgOfCell(ME),ME)=I

   !Part 5:
    IF(NE/=0)Then
	 NEdgOfCell(NE) = NEdgOfCell(NE) + 1
     InxEdgOfCell(NEdgOfCell(NE),NE)=-I
    EndIF

 End Do

  
!Part 6:
 Do I=1,NC

   !Part 7:
    Do J1=1,NEdgOfCell(I)
       E1 = InxEdgOfCell(J1,I)

      !Part 8:
       IF( E1>0 )Then
        P2_E1 = IDS(4,E1)
       Else
        P2_E1 = IDS(3,-E1)
	   EndIF

      !Part 9:
	   Do J2=J1+1,NEdgOfCell(I)
          E2 = InxEdgOfCell(J2,I)

         !Part 10:
          IF( E2>0 )Then
           P1_E2 = IDS(3,E2)
          Else
           P1_E2 = IDS(4,-E2)
          EndIF

         !Part 11:
          IF( P2_E1==P1_E2 )Then
	       E                 = InxEdgOfCell(J1+1,I)
	       InxEdgOfCell(J1+1,I) = InxEdgOfCell(J2,I)
           InxEdgOfCell(J2,I)   = E
          EndIF

       End Do

    End Do

 End Do 

!*********************************************************************************************
 End
!###########################################################################################











