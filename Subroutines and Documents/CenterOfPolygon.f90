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
Subroutine CenterOfPolygon(Dim,N,Points,X,Y,Cx,Cy) !--->>> Implementation of Method proposed by Paul Bourke (1988): 'Calculating The Area And Centroid Of A Polygon'
Implicit None
!===========================================================================================
Intent(In)::Dim,N,Points,X,Y
Intent(InOut)::Cx,Cy

Integer::Dim,N,Pi,Pi1,I
Integer,Dimension(1:N)::Points
Real(8)::Cx,Cy,A
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Notice: List of points comprising borders of the polygon must be consecutive(CW or CCW).
!--->>> NOTE: The geometric centroid of a convex object always lies in the object. <<<---

!Part 1:

A = 0.0

Cx = 0.0
Cy = 0.0

!Part 2:

do I=1,N
	
	Pi = Points(I)

	if(I+1 > N) then

		Pi1 = Points(1)

	else

		Pi1 = Points(I+1)

	endif

	A = A + X(Pi)*Y(Pi1) - X(Pi1)*Y(Pi)

	Cx = Cx + (X(Pi) + X(Pi1))*(X(Pi)*Y(Pi1) - X(Pi1)*Y(Pi))
	Cy = Cy + (Y(Pi) + Y(Pi1))*(X(Pi)*Y(Pi1) - X(Pi1)*Y(Pi))

end do

A = A/2

Cx = Cx/(6*A)
Cy = Cy/(6*A)

!===========================================================================================
End Subroutine CenterOfPolygon 
!*********************************************************************************************
