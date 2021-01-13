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
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function hasAngleLessThanSixDegrees(Dim,Corn,X,Y,Elm)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Elm

Integer::Dim,Elm,A,B,D,I,index
Integer,Dimension(1:Dim,1:4)::Corn
Logical::hasAngleLessThanSixDegrees
Real(8)::PI,Threshold,GetAngle
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
PI = 3.14159265358979d0

hasAngleLessThanSixDegrees = .False.

Threshold = PI/30 !-------------------------- 6 Degrees ------------------------------------
!Part 1:
do I=1,4

	 A = Corn(Elm,I)
	 
	 Call getNextCorner(Dim,Corn,Elm,A,B,index)
	 Call getPrevCorner(Dim,Corn,Elm,A,D,index)

	 if(GetAngle(Dim,A,B,D,X,Y) < Threshold) then

		hasAngleLessThanSixDegrees = .True.
		exit

	 endif

end do

!===========================================================================================
End Function hasAngleLessThanSixDegrees
!*********************************************************************************************
