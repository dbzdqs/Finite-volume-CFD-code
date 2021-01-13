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
 Subroutine Write_ConservativeVariables3D(Dim,NC,WNP1)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,NC,WNP1

 Integer::Dim,I,NC
 Real(8),Dimension(1:5,1:Dim)::WNP1
!**********************************************************************************************	
!Part 1:
 Open(1,File='ConservativeVariables.txt')
 Write(1,*) 'Ro  RU  RV  RW  RE '
 
 !Part 2:
 Do I=1,NC
	Write(1,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I),WNP1(5,I)
 End Do
 Close(1)
 
!**********************************************************************************************	
 END 
!##############################################################################################
 
 
 
 
