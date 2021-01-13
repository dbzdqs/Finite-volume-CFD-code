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
 Subroutine Write_CF(Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu)
 Implicit None
!**********************************************************************************************
 Intent(In)::Dim,Minf,Rinf,NFW1,NFW2,X,Y,IDS,DUY,Mu

 Integer::Dim,I,P1,P2,ME,NFW1,NFW2
 Real(8)::Minf,xm,ym,Rinf,CFC
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,DUY,Mu
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

!Part 5:
    Xm = 0.5*( X(P1)+X(P2) )

!Part 6:
    Write(21,*) Xm,  CFC *Mu(ME)*DUY(I)

 End do
 
 Close(21)
!**********************************************************************************************
 End
!##############################################################################################
 
 
 
 
 
