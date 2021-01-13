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
 Subroutine DeformableCell(Dim,Dead,Heir,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NDeformCell,IDeformCell,NP,NC)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dead,Heir,NPoint_Cell,IPoint_Cell,NConectCell,IConectCell,NP,NC
 Intent(Out  )::NDeformCell,IDeformCell
 
 Integer::Dim,I,J,Point,Cell,NDeformCell,Dead,Heir,NP,NC
 Integer,Dimension(1:Dim)::NPoint_Cell
 Integer,Dimension(1:Dim)::NConectCell
 Integer,Dimension(1:8,1:Dim)::IPoint_Cell
 Integer,Dimension(1:100,1:Dim)::IConectCell
 Integer,Dimension(1:100)::IDeformCell
!*********************************************************************************************
!Part 1:
 NDeformCell=0
 Do I=1,NConectCell(Dead) 
     Cell = IConectCell(I,Dead)
     
    !Part 2:
     Do J=1,NPoint_Cell(Cell)
         
        !Part 3:
		 Point = IPoint_Cell(J,Cell)
         IF( Point==Heir )Then
             NDeformCell = NDeformCell+1
             IDeformCell(NDeformCell) = Cell
             exit
         EndIF

     End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################
