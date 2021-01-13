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
!// Developed by: A. Hemati zadeh, Mechanical Eng., Amirkabir University of Technology     //!
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine BoundPointLabeling3D(Dim,Bound,NF,NP,IDS,FaceType)
 Implicit None
!*********************************************************************************************
 Intent(In   )::NF,NP,IDS,FaceType
 Intent(Out  )::Bound

 Integer::Dim,I,NF,NP
 Logical,Dimension(1:Dim)::Bound
 Integer,Dimension(1:6,1:Dim)::IDS
 Integer,Dimension(1:Dim)::FaceType
!*********************************************************************************************
!Part 1:
 Bound(:) = .FALSE.
 Do I=1,NF
     
   !Part 2:
    if( IDS(2,I)==0 )Then
     Bound(IDS(3,I))  = .TRUE.
     Bound(IDS(4,I))  = .TRUE.
     Bound(IDS(5,I))  = .TRUE.
     If(FaceType(I)==4)Then
      Bound(IDS(6,I)) = .TRUE.
     Endif   
    Endif
     
 EndDo
!*********************************************************************************************
 End
!###########################################################################################
