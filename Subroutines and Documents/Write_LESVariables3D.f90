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
!// Developed by: E. Akrami, Mechanical Eng., Amirkabir University of Technology           //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Write_LESVariables3D(Dim,NC,NP,Corn,X,Y,Z,time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean)
 Implicit None
!*********************************************************************************************
 Intent(In)::Dim,NC,NP,Corn,X,Y,Z,time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean

 Integer::Dim,J,NP,NC,Naverage
 Real(8)::time,umeanJ,vmeanJ,u2meanJ,v2meanJ,uvmeanJ
 Integer,Dimension(1:8,1:Dim)::Corn
 Real(8),Dimension(1:Dim)::X,Y,Z
 Real(8),Dimension(1:Dim)::sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean,uufluc,vvfluc,uvfluc
!*********************************************************************************************	
 Do J =1,NC

   !Part 1:
    umeanJ  =  sum_umean(J) /Naverage
    vmeanJ  =  sum_vmean(J) /Naverage 
    u2meanJ =  sum_u2mean(J)/Naverage 
    v2meanJ =  sum_v2mean(J)/Naverage
    uvmeanJ =  sum_uvmean(J)/Naverage

   !Part 2:
    uufluc(J) =  u2meanJ - umeanJ * umeanJ
    vvfluc(J) =  v2meanJ - vmeanJ * vmeanJ
    uvfluc(J) =  uvmeanJ - umeanJ * vmeanJ

 End Do   
     
!Part 3:
 Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,uufluc,"uufluc.plt")
 
 Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,vvfluc,"vvfluc.plt")
 
 Call Write_ScalarContour3D(Dim,NC,NP,X,Y,Z,Corn,uvfluc,"uvfluc.plt")

!Part 4:      
 Open(1,File='LESVariablesHistory.txt')
 Write(1,*) Naverage,time
 Do J=1,NC
    Write(1,*) sum_umean(J), sum_vmean(J), sum_u2mean(J), sum_v2mean(J), sum_uvmean(J)
 End Do
 close(1)

!**********************************************************************************************
 End
!##############################################################################################
