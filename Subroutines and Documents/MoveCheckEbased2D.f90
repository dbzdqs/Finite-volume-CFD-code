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
 Subroutine MoveCheckEbased2D(Dim,NF,NC,IDS,X,Y,Move)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NF,NC,IDS,X,Y
 Intent(Out  )::Move

 Integer::Dim,I,P1,P2,P3,NC,NF,ME,NE,Move
 Real(8)::SumX,SumY,DArea,DX,DY,DL
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::X,Y,A
!*********************************************************************************************	
!Part 1:
 Move=1

!Part 2:
 DO I=1,NC
    A(I)  = 0.0
 End do

!Part 3:
 DO I=1,NF

   !Part 4:
    ME = IDS(1,I)
    NE = IDS(2,I)
	P1 = IDS(3,I)
    P2 = IDS(4,I)

   !Part 5:
    DArea = X(P1)*Y(P2) - X(P2)*Y(P1)

   !Part 6:
    A(ME) = A(ME) + DArea / 2

   !Part 7:
    IF(NE/=0)A(NE) = A(NE) - DArea / 2

 End do

!Part 8:
 DO I=1,NC
    IF( A(I)<0. )Then
	 Move=-1
	 Exit
	EndIF
 End do
!*********************************************************************************************
 End
!###########################################################################################
