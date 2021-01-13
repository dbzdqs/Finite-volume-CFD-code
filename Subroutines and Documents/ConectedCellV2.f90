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
!// Date: May., 15, 2016                                                                   //!
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
 Subroutine ConectedCellV2(Dim,NC,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NP)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NPoint_Cell,IPoint_Cell,NP
 Intent(Out  )::NConectCell,IConectCell
 
 Integer::Dim,NC,I,J,Point,Temp,Cell,NP
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:Dim)::NPoint_Cell
 Integer,Dimension(1:4,1:Dim)::IPoint_Cell
 Integer,Dimension(1:50,1:Dim)::IConectCell
!*********************************************************************************************
!Part 1:
 NConectCell(1:NP) = 0

!Part 2:
 Do I=1,NC 
    Do J=1,NPoint_Cell(I)
        
      !Part 3:
       Point = IPoint_Cell(J,I)
       NConectCell(Point) = NConectCell(Point) + 1
	   IConectCell( NConectCell(Point),Point ) = I
       
    End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################