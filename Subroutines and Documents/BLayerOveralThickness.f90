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
!// Developed by: R. Amery, Mathmatical, Amirkabir university of Technology                //! 
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BLayerOveralThickness(Dim,OveralThick,StrmThikRatio,Xref,XBL,NPtCurv,NCurv,BLThick)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,OveralThick,StrmThikRatio,Xref,XBL,NPtCurv,NCurv
 Intent(Out  )::BLThick

 Integer::Dim,I,J,Sum
 Real(8)::Xdis
 Integer::NPtCurvs
 Integer::NCurv
 Real(8),Dimension(1:Dim)::XBL
 Real(8),Dimension(1:Dim)::BLThick
 Real(8),Dimension(1:100)::OveralThick
 Real(8),Dimension(1:100)::StrmThikRatio
 Real(8),Dimension(1:100)::Xref
 Integer,Dimension(1:100)::NPtCurv
!*********************************************************************************************
!Part 1:
 Sum = 0
 Do J=1,NCurv
    Do I=Sum+1,Sum+NPtCurv(J)
       BLThick(I) = OveralThick(J)
       XDis = dabs(Xref(J)-XBL(I))
       BLThick(I) = BLThick(I) + StrmThikRatio(J)*dsqrt(XDis)
    End Do
    Sum=Sum+NPtCurv(J)
 End Do
!*********************************************************************************************
 End
!###########################################################################################