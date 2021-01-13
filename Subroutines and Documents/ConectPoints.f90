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
!// Developed by: R. Rabii, Mechanical Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine ConectPoints(Dim,NF,IDS,NConectPoints,IConectPoints)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,IDS
 Intent(Out  )::NConectPoints,IConectPoints

 Integer::Dim,I,NF,P1,P2
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::NConectPoints
 Integer,Dimension(1:10,1:Dim)::IConectPoints
!*********************************************************************************************
!Part 1:
 Do I=1,NF
     
  !Part 2:
   P1 = IDS(3,I)
   P2 = IDS(4,I)
   
  !Part 3:
   NConectPoints(P1) = NConectPoints(P1) + 1
   NConectPoints(P2) = NConectPoints(P2) + 1
   
  !Part 4:
   IConectPoints( NConectPoints(P1), P1) = P2
   IConectPoints( NConectPoints(P2), P2) = P1

 End Do
!*********************************************************************************************
 End
!###########################################################################################
