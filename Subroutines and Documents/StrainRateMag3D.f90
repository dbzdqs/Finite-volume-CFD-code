
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
function StrainRateMag(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ)
 Implicit none

 Real(8)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ
 Real(8)::Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz
 Real(8)::SijSij,StrainRateMag

 Sxx = 0.5 * (DUX+DUX) ; Sxy = 0.5 * (DUY+DVX) ; Sxz = 0.5 * (DUZ+DWX)  
 Syx = Sxy             ; Syy = 0.5 * (DVY+DVY) ; Syz = 0.5 * (DVZ+DWY)
 Szx = Sxz             ; Szy = Syz             ; Szz = 0.5 * (DWZ+DWZ)

 SijSij  = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + &
           Syx*Syx + Syy*Syy + Syz*Syz + &
           Szx*Szx + Szy*Szy + Szz*Szz
    
 StrainRateMag = Dsqrt(2*SijSij)
    end 
    
 