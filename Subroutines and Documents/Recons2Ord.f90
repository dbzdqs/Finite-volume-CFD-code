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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Recons2Ord(Dim,J,I,Fc,GWNP1,Phi,DX,DY,Fc_face)
 Implicit None
!*********************************************************************************************	
 Intent(In  )::Dim,J,I,Fc,GWNP1,Phi,DX,DY
 Intent(Out )::Fc_face

 Integer::Dim,I,J
 Real(8)::DFc_Dx,DFc_Dy,Fc,Phi,DX,DY,Fc_face  
 Real(8),Dimension(1:2,1:4,1:Dim)::GWNP1
!*********************************************************************************************
!Part 1:
 DFc_Dx   = GWNP1(1,J,I)
 DFc_Dy   = GWNP1(2,J,I)

!Part 2:
 Fc_face = Fc + Phi * (DFc_Dx*DX + DFc_Dy*DY)   
!*********************************************************************************************    
 End 
!###########################################################################################