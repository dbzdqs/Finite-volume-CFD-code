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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function IsDiagonalInsideQuad(Dim,Corn,X,Y,E,A,C) !-- Implementation of 'Jordan Curve Theorem',Reference:http://erich.realtimerendering.com/ptinpoly/
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,E,A,C

Integer::Dim,E,A,B,C,Ai,Bi,Ci,Crossing,I
Integer,Dimension(1:Dim,1:4)::Corn
Logical::IsDiagonalInsideQuad
Real(8)::X_beg_ray,Y_beg_ray,X_end_ray,Y_end_ray,Dx0,Dy0,Dx1,Dy1,Deltax,Deltay,s,t,temp
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!----- Notice: A = N_c and C = N_d and AC is Diagonal to be Checked ------------------------
Crossing = 0
!Part 1:
do I=1,4
	if(Corn(E,I) == A) Ai = I
	if(Corn(E,I) == C) Ci = I	
end do

if(Ai == 1 .Or. Ai == 3) then
	
	if(Ai == 1) then
		Bi = 2
		B = Corn(E,2)
	else
		Bi = 4
		B = Corn(E,4)
	endif

elseif(Ai == 2 .Or. Ai == 4) then
	
	if(Ai == 2) then
		Bi = 3
		B = Corn(E,3)
	else
		Bi = 1
		B = Corn(E,1)
	endif

endif
!Part 2:
!------------- Coordinate of midpoint of the diagonal as beginning of RAY ------------------

X_beg_ray = (X(A) + X(C))/2
y_beg_ray = (Y(A) + Y(C))/2

!------------ Coordinate of midpoint of an optional edge as end of RAY ----------------------

X_end_ray = (X(A) + X(B))/2
y_end_ray = (Y(A) + Y(B))/2

!Part 3:

!---------------------------------- RAY in Vector form -------------------------------------

Dx0 = X_end_ray - X_beg_ray
Dy0 = Y_end_ray - Y_beg_ray

do I=1,4
	
	if(I < 4) then

		Dx1 = X(Corn(E,I+1)) - X(Corn(E,I))
		Dy1 = Y(Corn(E,I+1)) - Y(Corn(E,I))

		Deltax = X(Corn(E,I)) - X_beg_ray
		Deltay = Y(Corn(E,I)) - Y_beg_ray 

		temp = Dx0*Dy1 - Dx1*Dy0
		
		if(temp /= 0) then

			s = (Deltax*Dy1 - Deltay*Dx1)/temp !--- s Corresponds to Ray that we suppose it an INFINITE RAY ---
			t =	(Deltax*Dy0 - Deltay*Dx0)/temp

			if(s>=0 .And. t>=0 .And. t<=1) then
				Crossing = Crossing + 1
			endif
		
		endif 
		  
	else

		Dx1 = X(Corn(E,1)) - X(Corn(E,4))
		Dy1 = Y(Corn(E,1)) - Y(Corn(E,4))

		Deltax = X(Corn(E,4)) - X_beg_ray
		Deltay = Y(Corn(E,4)) - Y_beg_ray 

		temp = Dx0*Dy1 - Dx1*Dy0
		
		if(temp /= 0) then

			s = (Deltax*Dy1 - Deltay*Dx1)/temp !--- s Corresponds to Ray that we suppose it an INFINITE RAY ---
			t =	(Deltax*Dy0 - Deltay*Dx0)/temp

			if(s>=0 .And. t>=0 .And. t<=1) then
				Crossing = Crossing + 1
			endif
		
		endif

	endif

end do
!Part 4:
if(MOD(Crossing,2) /= 0) then
	IsDiagonalInsideQuad = .True.
else
	IsDiagonalInsideQuad = .False.
endif 

!===========================================================================================
End Function IsDiagonalInsideQuad
!*********************************************************************************************
