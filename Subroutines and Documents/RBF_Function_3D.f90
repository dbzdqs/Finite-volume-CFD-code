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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Function3D(DL,Temp)
 Implicit None
!********************************************************************************************* 
 Intent(In  )::DL
 Intent(Out )::Temp

 Real(8)::Temp,DL
!********************************************************************************************* 

!Temp=((1-(DL/5.))**4)*((4*(DL/5.))+1)   !toop

!Temp=Dl
!Temp=(1-(DL/2.))**2                      !wing
 
!Temp=((1-DL/2.)**6)*((35*((Dl/2.)**2))+(18*Dl/2.)+3)
!Temp=((1-DL)**8)*((32*(Dl**3))+(25*(Dl**2))+(8*Dl)+1)
!Part 1:

!IF(DL>1.0e-8)Then    
! Temp = (DL*DL)*(dlog(DL))
!Else 
!Temp = 0.0       
!End If    
 
!Part 2:
!Temp =1./( 1. + (DL/2.)**2)**2              !wing bending-NSTP=40-H=50   **
Temp =1./( 1. + (DL)**2)             !wing bending-NSTP=40-H=50   **

!Part 3:
!Temp = 1. / (1. + (DL/1.)**2)**2            ! NACA0012 Plunging 

!Temp = 1. / (1. + (DL/6.)**2)**2           !   beam bending -NSTP=30-H=100

!Temp = 1. / (1. + (DL/4.)**2)**2           !   ball plunging -NSTP=40-H=10  **

!Temp = 1. / (1. + (DL/5.)**2)**2           !   wing bending -NSTP=30-H=50

!Temp = 1. / (1. + (DL/2.)**2)**2           !   wing plunging1 -NSTP=40-H=15  **
!Temp = 1. / (1. + (DL/4.)**2)**2           !   wing plunging2 -NSTP=40-H=15  **

!Part 4:
!Temp = exp( -DL*DL )                           !   wing bending -NSTP=30-H=50
!*********************************************************************************************
 End
!###########################################################################################
 