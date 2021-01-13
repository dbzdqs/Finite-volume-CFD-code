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
!// Developed by: M. Valadkhani, Mechanical Eng., Amirkabir University of Technology       //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BC_Wall_ALE(Dim,NFW1,NFW2,IDS,GM,WNP1,P,Face_Velocity,WB)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NFW1,NFW2,IDS,GM,WNP1,P,Face_Velocity
 Intent(Out  )::WB

 Integer::Dim,I,NFW1,NFW2,ME
 Real(8)::GM,GM1,PB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::P
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
  Real(8),Dimension(1:2,1:Dim)::Face_Velocity
!*********************************************************************************************	
!Part 1:  
 GM1= GM-1.
 
!Part 2:
 DO I=NFW1+1,NFW2
     
    !Part 3:
     ME = IDS(1,I)
     
    !Part 4:
     PB = P(ME)
     
    !Part 5:
     WB(1,I) = WNP1(1,ME)
     
    !Part 6:
    !Rezaei's Comment
    !ALE Majid 
    !WB(2,I) = 0.0 + WB(1,I)*Face_Velocity(1,I)
    !WB(3,I) = 0.0 + WB(1,I)*Face_Velocity(2,I)
     WB(2,I) = Face_Velocity(1,I)
     WB(3,I) = Face_Velocity(2,I)
     
    !Part 7:
     WB(4,I) = PB/GM1 + 0.5*( WB(2,I)* WB(2,I) + WB(3,I)*WB(3,I))/WB(1,I)
     WB(5,I) = PB
     
 End do
!*********************************************************************************************
 End
!###########################################################################################
