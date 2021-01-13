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
!// Date: Mar., 05, 2013                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CalculateConstants(Dim,NBE,BFP,X,Y,DELTA,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,X,Y
Intent(Out)::DELTA,COINCIDENT_TOLERANCE

Integer::Dim,NBE,I
Integer,Dimension(1:Dim,1:2)::BFP
Real(8)::DELTA,COINCIDENT_TOLERANCE,L,SE,LE,GetNorm
Real(8),Dimension(1:Dim)::X,Y
!=========================================================================================== 

!-------------- Finding Smallest Edge (SE) and Longest Edge (LE) on the boundary -----------

!Part 1:

SE = GetNorm(Dim,BFP(1,1),BFP(1,2),X,Y)
LE = SE

DO I=1,NBE
 
	L = GetNorm(Dim,BFP(I,1),BFP(I,2),X,Y)

	IF(LE < L) LE = L
	IF(SE > L) SE = L


END DO

!Part 2:

COINCIDENT_TOLERANCE = SE*0.0001

DELTA = LE*0.00001

!===========================================================================================
End Subroutine CalculateConstants
!*********************************************************************************************
