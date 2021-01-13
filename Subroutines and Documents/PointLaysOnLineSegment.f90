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
Function PointLaysOnLineSegment(Dim,X,Y,A,B,P) ! -->> Determines whether a point is on a LINE SEGMENT or not <<--
Implicit None
!===========================================================================================
Intent(In)::Dim,X,Y,A,B,P

Integer::Dim,A,B,P
Logical::PointLaysOnLineSegment
Real(8)::dx1,dy1,dx2,dy2,cross
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: https://stackoverflow.com/questions/11907947/how-to-check-if-a-point-lies-on-a-line-between-2-other-points: answered by: AnT 
!Part 1:
!------------ Vector from A to P ---------
dx1 = X(P) - X(A)
dy1 = Y(P) - Y(A)

!------------ Vector from A to B ---------
dx2 = X(B) - X(A)
dy2 = Y(B) - Y(A)

!Part 2:

cross = dx1*dy2 - dx2*dy1 

if(cross /= 0) then

	PointLaysOnLineSegment = .False.

else

!Part 3:    
    
	if(ABS(dx2) >= ABS(dy2)) then

		if(dx2 > 0) then

			if(X(A) <= X(P) .And. X(P) <= X(B)) then

				PointLaysOnLineSegment = .True.

			else

				PointLaysOnLineSegment = .False.

			endif

		else

			if(X(B) <= X(P) .And. X(P) <= X(A)) then

				PointLaysOnLineSegment = .True.

			else

				PointLaysOnLineSegment = .False.

			endif

		endif

	else

		if(dy2 > 0) then

			if(Y(A) <= Y(P) .And. Y(P) <= Y(B)) then

				PointLaysOnLineSegment = .True.

			else
				
				PointLaysOnLineSegment = .False.

			endif

		else

			if(Y(B) <= Y(P) .And. Y(P) <= Y(A)) then

				PointLaysOnLineSegment = .True.

			else

				PointLaysOnLineSegment = .False.

			endif

		endif

	endif

endif

!===========================================================================================
End Function PointLaysOnLineSegment
!*********************************************************************************************
