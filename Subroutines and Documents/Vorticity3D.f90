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
 function Vorticity(DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ)
 Implicit none

 Real(8)::DUX,DUY,DUZ,DVX,DVY,DVZ,DWX,DWY,DWZ
 Real(8)::Omgaxy,Omgaxz,Omgayx,Omgayz,Omgazx,Omgazy
 Real(8)::OmgaijOmgaij,Vorticity

                            Omgaxy = 0.5 * (DUY-DVX) ; Omgaxz = 0.5 * (DUZ-DWX)  
 Omgayx = 0.5 * (DVX-DUY)                            ; Omgayz = 0.5 * (DVZ-DWY)
 Omgazx = 0.5 * (DWX-DUZ) ; Omgazy = 0.5 * (DWY-DVY)            

 OmgaijOmgaij  =                 Omgaxy*Omgaxy + Omgaxz*Omgaxz + &
                 Omgayx*Omgayx +                 Omgayz*Omgayz + &
                 Omgazx*Omgazx + Omgazy*Omgazy
    
 Vorticity = Dsqrt(2*OmgaijOmgaij)
 end 