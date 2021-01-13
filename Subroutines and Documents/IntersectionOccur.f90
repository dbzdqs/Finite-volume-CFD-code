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
Function IntersectionOccur(Dim,Nc,Nd,c1,c2,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,Nc,Nd,c1,c2,X,Y

Integer::Dim,Nc,Nd,c1,c2
Logical::IntersectionOccur
Real(8)::Dx0,Dx1,Dy0,Dy1,Deltax,Deltay,s,t,temp
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
if((c1==Nc .Or. c1==Nd) .Or. (c2==Nc .Or. c2==Nd)) then
	IntersectionOccur = .False.
else
!Part 2:
	Dx0 = X(Nd) - X(Nc)
	Dy0 = Y(Nd) - Y(Nc)

	Dx1 = X(c2) - X(c1)
	Dy1 = Y(c2) - Y(c1)
!Part 3:
	Deltax = X(c1) - X(Nc)
	Deltay = Y(c1) - Y(Nc)

	temp = Dx0*Dy1 - Dx1*Dy0

	if(temp /= 0) then

		s = Real((Deltax*Dy1 - Dx1*Deltay)/temp,8)
		t = Real((Deltax*Dy0 - Dx0*Deltay)/temp,8)

		if(s>0 .And. s<1 .And. t>0 .And. t<1) then
			
			IntersectionOccur = .True.

		else

			IntersectionOccur = .False.

		endif

	else

		IntersectionOccur = .False.

	endif

endif
!===========================================================================================
End Function IntersectionOccur
!*********************************************************************************************
