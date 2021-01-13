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
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Read_CST_Parameters(BSOrder_Up,BSOrder_Lw,UpRegion,LwRegion)
 Implicit None
!*********************************************************************************************
 Intent(Out  )::BSOrder_Up,BSOrder_Lw,UpRegion,LwRegion

 Integer::BSOrder_Up,BSOrder_Lw,UpRegion,LwRegion
!*********************************************************************************************
!Part 1:
 Open(40, FILE='CST_Parameters.txt', STATUS='Old')

 Read(40,*) UpRegion , BSOrder_Up
 Read(40,*) LwRegion , BSOrder_Lw

 Close(40)
!*********************************************************************************************
 End
!###########################################################################################