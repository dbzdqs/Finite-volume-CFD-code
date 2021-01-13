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
 Subroutine ResCFL_Calculator(Dim,NC,Rn,Ncyc,First_Res,CFL,Res)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,Rn,Ncyc,First_Res
 Intent(InOut)::CFL,Res

 Integer::Dim,Ncyc,NC,I,J
 Real(8)::Res,First_Res,CFL,time_end
 Real(8),Dimension(1:4,1:Dim)::Rn    
!*********************************************************************************************	
 Open(111,File='ResMass_Iter.Plt')
 Open(110,File='ResMass_Time.Plt')

!Part 1:
 Res = 0.0
 Do I=1,NC
    Do J=1,4
       Res = Res + Rn(J,I)*Rn(J,I)
    End Do
 End Do    
 Res = Sqrt(Res)
  
!Part 2:  
 CFL = 10.0 - 22.0*Log(Res/First_Res)     

!Part 3:
 Call CPU_Time(time_end)
 Write(110,*) time_end,Res
 Write(111,*) Ncyc,Res
!*********************************************************************************************
 End
!###########################################################################################