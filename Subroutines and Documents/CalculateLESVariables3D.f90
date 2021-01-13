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
!// Date: May., 15, 2016                                                                   //!
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
 Subroutine CalculateLESVariables3D(Dim,Init,Ncyc,NC,WNP1,time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,WNP1,Init,Ncyc
 Intent(InOut)::time,Naverage,sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean
 
 Integer::Dim,J,NC,Naverage,Init,Ncyc
 Real(8)::U,V,time
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::sum_umean,sum_vmean,sum_uvmean,sum_u2mean,sum_v2mean
!*********************************************************************************************	
!Part 1:
 IF(Init==0 .and. Ncyc==1)Then
    
  Do J =1,NC
     sum_umean(J)  = 0.0
     sum_vmean(J)  = 0.0
     sum_u2mean(J) = 0.0
     sum_v2mean(J) = 0.0
     sum_uvmean(J) = 0.0
  End Do  

!Part 2:
 ElseIF(Init==1 .and. Ncyc==1)Then
    
  Open(1,File='LESVariablesHistory.txt')
  Read(1,*) Naverage,time
  Do J=1,NC
     Read(1,*) sum_umean(J), sum_vmean(J), sum_u2mean(J), sum_v2mean(J), sum_uvmean(J)
  End Do
  close(1)

 EndIF
 
!Part 3:
 Do J =1,NC
    U = WNP1(2,J) / WNP1(1,J)
    V = WNP1(3,J) / WNP1(1,J)
        
    sum_umean(J) =  sum_umean(J) + U 
    sum_vmean(J) =  sum_vmean(J) + V 
        
    sum_u2mean(J) = sum_u2mean(J) + U*U 
    sum_v2mean(J) = sum_v2mean(J) + V*V
        
    sum_uvmean(J) = sum_uvmean(J) + U*V
 End Do     
  
!Part 4:
 Naverage = Naverage + 1
!**********************************************************************************************
 End
!##############################################################################################
