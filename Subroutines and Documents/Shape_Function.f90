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
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Shape_Function(DimU,DimL,NPtCurv_Up,X_Up,Y_Up,NPtCurv_Lw,X_Lw,Y_Lw,Zita_TE,S_Up,S_Lw)
Implicit None
!********************************************************************************************* 
 Intent(In   )::DimU,DimL,NPtCurv_Up,X_Up,Y_Up,NPtCurv_Lw,X_Lw,Y_Lw,Zita_TE
 Intent(Out  )::S_Up,S_Lw

 Integer::DimU,DimL,I,J,NPtCurv_Up,NPtCurv_Lw
 Real(8)::Zita_TE 
 Real(8),Dimension(1:DimU)::X_up,Y_up,S_up
 Real(8),Dimension(1:DimL)::X_Lw,Y_Lw,S_Lw    
!*********************************************************************************************
!Part1:
 Do I=2,NPtCurv_Up-1
    S_up(I) = (Y_up(I) - X_up(I)*Zita_TE)/(SQRT(X_up(I)) * (1. - X_up(I)))
 End Do

!Part2:
 Do J=2,NPtCurv_Lw-1
    S_Lw(J) = (Y_Lw(J) + X_Lw(J)*Zita_te)/(SQRT(X_Lw(J)) * (1. - X_Lw(J)))
 End Do
!*********************************************************************************************
 End
!###########################################################################################