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
 Subroutine Kw_Init(Dim,NC,MR,WTNP1,Mut) 
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,MR
 Intent(Out  )::WTNP1,Mut

 Integer::Dim,J,NC
 Real(8)::MR
 Real(8),Dimension(1:2,1:Dim)::WTNP1
 Real(8),Dimension(1:Dim)::Mut
!*********************************************************************************************	
!Part 1: 	
 Do J=1,NC
    WTNP1(1,J)=1.0e-6
    WTNP1(2,J)=5.0
    Mut(J) = (1/MR)*(WTNP1(1,J)/WTNP1(2,J)) 
 End Do

!*********************************************************************************************
 End
!###########################################################################################
