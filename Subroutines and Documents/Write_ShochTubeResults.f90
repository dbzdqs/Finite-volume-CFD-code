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
 Subroutine Write_ShockTubeResults(Dim,NC,Xc,WNP1,P)
 Implicit None
!*********************************************************************************************
 Intent(In)::Dim,NC,Xc,WNP1,P

 Integer::Dim,I,NC
 Real(8),Dimension(1:Dim)::P,Xc,Yc
 Real(8),Dimension(1:4,1:Dim)::WNP1
!*********************************************************************************************	
!Part 1:

 Open(20,File='ShockTubeResults.Plt')


 WRITE(20,*)'VARIABLES="X","Density","Velocity","Energy","Press"'  !
 WRITE(20,*)'ZONE'
!Part 4:
 Do I=1,NC
	Write(20,*) Xc(I)  ,WNP1(1,I) ,WNP1(2,I)/WNP1(1,I),WNP1(4,I)/WNP1(1,I),P(I)
 End Do

 Close(20)
!**********************************************************************************************
 End
!##############################################################################################
