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
 Subroutine	GCL(Dim,NF,IDS,GF,DTmin,A)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NF,IDS,GF,DTmin
 Intent (InOut  )::A

 Integer::Dim,I,NF,ME,NE
 Real(8)::DTmin
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::A,GF
!*********************************************************************************************
!Part 1:
 !DO I=NF1+1,NF2
 !Majid
 DO I=1,NF
     
    !Part 2:
     ME = IDS(1,I)
     NE = IDS(2,I)
     
    !Part 3:
     A(ME) = A(ME) + DTmin*GF(I)
     
    !Part 4:
    If(NE .NE. 0)  A(NE) = A(NE) - DTmin*GF(I)
     
 End do

!*********************************************************************************************
 End
!########################################################################################### 
