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
 Subroutine Write_CP(Dim,GM,Minf,NFW1,NFW2,X,Y,IDS,P)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,GM,Minf,NFW1,NFW2,X,IDS,P

 Integer::Dim,I,P1,P2,ME,NFW1,NFW2
 Real(8)::Xm,GM,Minf
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,P
!**********************************************************************************************	
!Part 1:
 Open(2,File='CP.Plt')

!Part 2:
 DO I=NFW1+1,NFW2

   !Part 3:
    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)

    !Part 4:
    Xm = 0.5*( X(P1)+X(P2))

   !Part 5:
	Write(2,*) Xm,(GM*P(ME)-1.0)/(0.5*GM*Minf*Minf)

 End do
 
 Close(2)
!**********************************************************************************************
 End
!##############################################################################################
 
 
 
 
 
