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
 Subroutine Write_CF3D(Dim,Minf,Rinf,NFW1,NFW2,X,Y,Z,IDS,DUY,Mu)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,Minf,Rinf,NFW1,NFW2,X,Y,Z,IDS,DUY,Mu

 Integer::Dim,I,P1,P2,P3,P4,ME,NFW1,NFW2
 Real(8)::Minf,xm,ym,zm,Rinf,CFC
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,Z,DUY,Mu
!**********************************************************************************************	
!Part 1:
 Open(21,File='CF.Plt')


!Part 2: 
 CFC = 2/(Rinf*Minf)
 
 !Part 3:
 DO I=NFW1+1,NFW2
 
 !Part 4:

    ME = IDS(1,I)
    P1 = IDS(3,I)
    P2 = IDS(4,I)
    P3 = IDS(5,I)
    P4 = IDS(6,I)

!Part 5:
    Zm = ( Z(P1)+Z(P2)+Z(P3)+Z(P4) ) / 4.
    Xm = ( X(P1)+X(P2)+X(P3)+X(P4) ) / 4.
    Ym = ( Y(P1)+Y(P2)+Y(P3)+Y(P4) ) / 4.

!Part 6:
    IF( Zm>0.0 .and. Zm<0.6 ) Write(21,*) Xm,  CFC *Mu(ME)*DUY(I)

 End do
 
 Close(21)
!**********************************************************************************************
 End
!##############################################################################################
 
 
 
 
 
