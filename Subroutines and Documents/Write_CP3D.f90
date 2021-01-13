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
 Subroutine Write_CP3D(Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,Minf,NFW1,NFW2,X,Y,Z,IDS,P

 Integer::Dim,I,P1,P2,P3,P4,ME,NFW1,NFW2
 Real(8)::Minf,xm,ym,zm
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,Z,P
!**********************************************************************************************	
!Part 1:
 Open(2,File='CP.Plt')

!Part 2:
 DO I=NFW1+1,NFW2


   !Part 3:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
    P3 = IDS(5,I)
    P4 = IDS(6,I)

    !Part 4:
    Zm = ( Z(P1)+Z(P2)+Z(P3)+Z(P4) ) / 4.
    Xm = ( X(P1)+X(P2)+X(P3)+X(P4) ) / 4.
    Ym = ( Y(P1)+Y(P2)+Y(P3)+Y(P4) ) / 4.

   !Part 5:
    IF( Zm>0.0 .and. Zm<3.6 )Then
	Write(2,*) Xm,(1.4*P(ME)-1.0)/(0.5*1.4*Minf*Minf)
	Endif

 End do
 
 Close(2)
!**********************************************************************************************
 End
!##############################################################################################
 
 
 
 
 
