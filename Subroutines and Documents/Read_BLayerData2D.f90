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
 Subroutine Read_BLayerData2D(Dim,NR,Region_BL,Dis1,OveralThick,Ratio,StrmThikRatio,Xref,NBL,DistributionType)
 Implicit None
!*********************************************************************************************
 Intent (In   )::Dim,NR
 Intent (Out  )::Region_BL,Dis1,OveralThick,Ratio,StrmThikRatio,Xref,NBL,DistributionType

 Integer::Dim,I,NR,NBL,J
 Integer,Dimension(1:100)::Region_BL
 Integer,Dimension(1:100)::DistributionType
 Real(8),Dimension(1:100)::OveralThick,Dis1,Ratio,StrmThikRatio,Xref
!*********************************************************************************************
!Part 1:
 Open(1,File='BLData2D.Txt')

!Part 2:
 DO I=1,NR
    Read(1,*) Region_BL(I)
 END DO

!Part 3:
 Read(1,*) NBL
 
!Part 4:
 Read(1,*)
 Read(1,*)
 DO I=1,NR
     Read(1,*) Xref(I),OveralThick(I), Dis1(I),Ratio(I),StrmThikRatio(I),DistributionType(I)
 END DO
 
 Close(1)
!*********************************************************************************************
 END 
!###########################################################################################