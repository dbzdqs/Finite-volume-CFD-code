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
 Subroutine ResMass3D(Dim,NC,WNP1,WN,DT,Ncyc,Rm)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,WNP1,WN,Ncyc,DT
 Intent(Out  )::Rm

 Integer::Dim,I,Ncyc,NC
 Real(8)::Rm,time_end
 Real(8),Dimension(1:5,1:Dim)::Wnp1,WN
 Real(8),Dimension(1:Dim)::DT
!*********************************************************************************************	
!Part 1:
 Open(111,File='ResMass_Iter.Plt')
 Open(110,File='ResMass_Time.Plt')
 	
!Part 2:
 Rm = 0.0

!Part 3:
 DO I=1,NC
    Rm = Rm + Abs( Wnp1(1,I)-WN(1,I) ) / DT(i)
 End Do

!Part 4:
 Rm = Rm / NC

!Part 5:
 Rm = Dlog10(Rm+1.0E-16)
 Write(111,*) Ncyc,Rm
 
 Call CPU_Time(time_end)
 Write(110,*) time_end,Rm

!*********************************************************************************************
 End
!###########################################################################################