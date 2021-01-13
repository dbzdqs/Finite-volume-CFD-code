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
 Subroutine Read_BDisplacement(Dim,NStp,NBP,IBP,DelX,DelY)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim
 Intent(Out  )::NStp,NBP,IBP,DelX,DelY

 Integer::Dim,I,NBP,NStp
 Integer,Dimension(1:Dim)::IBP
 Real(8),Dimension(1:Dim)::DelX,DelY
!*********************************************************************************************
!Part 1:
 Open(11,File='BDisplacement.Txt')

!Part 2:
 Read(11,*) NStp

!Part 3:
 Read(11,*)
 Read(11,*) NBP

!Part 4:	  
 Do I=1,NBP
    Read(11,'(I15)',Advance='No') IBP(I) 
    Read(11,*)   DelX( IBP(I) ),DelY( IBP(I) )
 End Do
!*********************************************************************************************
 End
!###########################################################################################

