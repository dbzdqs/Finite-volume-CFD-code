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
!// Date: Apr., 25, 2013                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RotatePt(Xo,Yo,Tet,Xp,Yp)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Xo,Yo,Tet
 Intent(InOut)::Xp,Yp

 Real(8)::Xo,Yo,Xp,Yp,PI,Tet,Teta,X_Vec,Y_Vec,X_Rot,Y_Rot
!*********************************************************************************************
 
 PI = 4*Atan(1.)
 Teta = Tet * (PI/180.)
   
 X_Vec = Xp - Xo
 Y_Vec = Yp - Yo

 X_Rot = Xo + X_Vec*cos(Teta) - Y_Vec*sin(Teta)
 Y_Rot = Yo + X_Vec*sin(Teta) + Y_Vec*cos(Teta)
    
 Xp= X_Rot 
 Yp= Y_Rot 

!*********************************************************************************************
 End
!###########################################################################################