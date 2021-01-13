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
 Subroutine Velocity_CellGrad(Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB,DUX_C,DUY_C,DVX_C,DVY_C)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,NF1,NF2,IDS,A,NX,NY,WNP1,WB
 Intent(Out  )::DUX_C,DUY_C,DVX_C,DVY_C

 Integer::Dim,I,NC,NF,ME,NE,NF1,NF2
 Real(8)::U,V,Rho,Area,NXX,NYY
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::DUX_C,DUY_C,DVX_C,DVY_C,NX,NY,A
!********************************************************************************************* 
!Part 1:
 DO I=1,NC
    DUX_C(I) = 0.0
    DUY_C(I) = 0.0
    DVX_C(I) = 0.0
    DVY_C(I) = 0.0
 End Do

!Part 2:
 DO I=NF1+1,NF2

   !Part 3:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 4:
    U  = 0.5*( WNP1(2,ME)/WNP1(1,ME)+WNP1(2,NE)/WNP1(1,NE) )
    V  = 0.5*( WNP1(3,ME)/WNP1(1,ME)+WNP1(3,NE)/WNP1(1,NE) )

   !Part 5:
    NXX = NX(I)
    NYY = NY(I)

    DUY_C(ME)  = DUY_C(ME)  + U*NYY
    DUX_C(ME)  = DUX_C(ME)  + U*NXX
    DVX_C(ME)  = DVX_C(ME)  + V*NXX
    DVY_C(ME)  = DVY_C(ME)  + V*NYY

    DUY_C(NE)  = DUY_C(NE)  - U*NYY
    DVY_C(NE)  = DVY_C(NE)  - V*NYY
    DUX_C(NE)  = DUX_C(NE)  - U*NXX
    DVX_C(NE)  = DVX_C(NE)  - V*NXX

 End Do

!Part 6:
 DO I=NF2+1,NF

   !Part 7:
    ME = IDS(1,I)

    NXX = NX(I)
    NYY = NY(I)

   !Part 8:
    Rho = WB(1,I)
    U   = WB(2,I) /Rho
    V   = WB(3,I) /Rho

   !Part 9:
    DUY_C(ME)  = DUY_C(ME)    + U*NYY
    DVY_C(ME)  = DVY_C(ME)    + V*NYY
    
    DUX_C(ME)  = DUX_C(ME)    + U*NXX
    DVX_C(ME)  = DVX_C(ME)    + V*NXX
 End Do 
 
!Part 10:
 DO I=1,NC
    Area = A(I)

    DUY_C(I) = DUY_C(I)/Area
    DUX_C(I) = DUX_C(I)/Area
    DVY_C(I) = DVY_C(I)/Area
    DVX_C(I) = DVX_C(I)/Area

               

  
 End Do
!*********************************************************************************************
 End
!###########################################################################################

