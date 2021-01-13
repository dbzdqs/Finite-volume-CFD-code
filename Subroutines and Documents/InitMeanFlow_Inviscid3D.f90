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
!// Date: April, 01, 2017                                                                  //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine InitMeanFlow_Inviscid3D(Dim,Init,NC,ALF,Minf,GM,R0,P0,C0,U0,V0,W0,WNP1)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Init,NC,Minf,ALF
 Intent(Out  )::GM,R0,P0,C0,U0,V0,W0,WNP1

 Integer::Dim,Init,I,NC
 Real(8)::PI,ALF,Minf,GM,R0,P0,T0,Ts,Tt,C0,U0,V0,W0,E0,ALFA,B0,U,V
 Real(8),Dimension(1:5,1:Dim)::WNP1
!*********************************************************************************************	
!Part 1:
 PI  = 4.0*Atan(1.0)
 ALFA= ALF*PI/180.

!Part 2:
 GM  = 1.4

!Part 3:
 R0 = 1.0

!Part 4:
 P0 = 1.0/GM  

!Part 5:
 T0 = GM*P0/R0  

!Part 6:
 C0 = SQRT(GM*P0/R0)

!Part 7:
 U0 = Minf*C0*COS(ALFA)
 W0 = 0.0
 V0 = Minf*C0*SIN(ALFA)

!Part 8:
 E0 = P0/(R0*(GM-1))+ 0.5*(U0*U0 + V0*V0 + W0*W0) 

!Part 9:
 DO I=1,NC
    WNP1(1,I) = R0
    WNP1(2,I) = R0*U0
    WNP1(3,I) = R0*V0
    WNP1(4,I) = R0*W0
    WNP1(5,I) = R0*E0
 end do

!Part 10:
 IF(Init==1)Then
  Open(1,File='ConservativeVariables.txt')
	 Read(1,*)
  Do I=1,NC
	 Read(1,*) WNP1(1,I),WNP1(2,I),WNP1(3,I),WNP1(4,I),WNP1(5,I)  !
  End Do       
  close(1)
 Endif
!*********************************************************************************************
 End
!###########################################################################################