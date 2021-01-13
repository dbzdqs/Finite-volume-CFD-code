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
Function SwapIsSafe(Dim,Corn,X,Y,tri_down,tri_up)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,tri_down,tri_up

Integer::Dim,tri_down,tri_up,P,Q,R,S,I,Iface1
Integer,Dimension(1:Dim,1:4)::Corn
Logical::SwapIsSafe
Real(8)::Dx0,Dy0,Dx1,Dy1,Deltax,Deltay,temp,t,u
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
if(Corn(tri_down,4) == 0 .And. Corn(tri_up,4) == 0) then
!Part 1:	
	do I=1,3
		if(Corn(tri_down,I) /= Corn(tri_up,1) .And. Corn(tri_down,I) /= Corn(tri_up,2) .And. Corn(tri_down,I) /= Corn(tri_up,3)) then
			if(I == 1) then
				P = Corn(tri_down,2)
				Q = Corn(tri_down,3)
			elseif(I == 2) then
				P = Corn(tri_down,3)
				Q = Corn(tri_down,1)
			elseif(I == 3) then
				P = Corn(tri_down,1)
				Q = Corn(tri_down,2)
			endif	
		endif
	end do
!Part 2:
	do I=1,3 
		if(Corn(tri_down,I) /= P .And. Corn(tri_down,I) /= Q) then
			Iface1 = I
			R = Corn(tri_down,I)
			exit
		endif
	end do

	do I=1,3 
		if(Corn(tri_up,I) /= P .And. Corn(tri_up,I) /= Q) then
			S = Corn(tri_up,I)
			exit
		endif
	end do
!Part 3:
	!-------------------------- D0: Vector from P to Q -----------------------------
	Dx0 = X(Q) - X(P)
	Dy0 = Y(Q) - Y(P)

	!-------------------------- D1: Vector from R to S -----------------------------
	Dx1 = X(S) - X(R)
	Dy1 = Y(S) - Y(R)
!Part 4:
	Deltax = X(R) - X(P)
	Deltay = Y(R) - Y(P)

	temp = Dx0*Dy1 - Dx1*Dy0 

	t = (Deltax*Dy1 - Deltay*Dx1)/temp 
	u = (Deltax*Dy0 - Deltay*Dx0)/temp

	if(t>0 .And. t<1 .And. u>0 .And. u<1) then
		SwapIsSafe = .True.
	else
		SwapIsSafe = .False.
	endif

else
	SwapIsSafe = .False.
endif
!===========================================================================================
End Function SwapIsSafe
!*********************************************************************************************
