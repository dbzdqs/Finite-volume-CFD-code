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
Function IsBoundaryElement(Dim,Corn,NBE,BFP,E)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,NBE,BFP,E

Integer::Dim,NBE,I,J,A,B,C,D,E
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn
Logical::IsOnExteriorBoundary,IsBoundaryElement
!===========================================================================================
do I=1,2
!Part 1:
	J = I + 2
	A = Corn(E,I)
	C = Corn(E,J)

	if(I == 1) then

		B = Corn(E,I + 1)
		D = Corn(E,J + 1)
		 
	elseif(I == 2) then

		B = Corn(E,I - 1)
		D = Corn(E,J - 1)

	endif
!Part 2:
	if(IsOnExteriorBoundary(Dim,NBE,BFP,A) .And. IsOnExteriorBoundary(Dim,NBE,BFP,B)) then

		IsBoundaryElement = .True.

	elseif(IsOnExteriorBoundary(Dim,NBE,BFP,A) .And. IsOnExteriorBoundary(Dim,NBE,BFP,D)) then

		IsBoundaryElement = .True.

	elseif(IsOnExteriorBoundary(Dim,NBE,BFP,B) .And. IsOnExteriorBoundary(Dim,NBE,BFP,C)) then

		IsBoundaryElement = .True.

	elseif(IsOnExteriorBoundary(Dim,NBE,BFP,C) .And. IsOnExteriorBoundary(Dim,NBE,BFP,D)) then

		IsBoundaryElement = .True.

	else

		IsBoundaryElement = .False.

	endif

end do
!===========================================================================================
End Function IsBoundaryElement 
!*********************************************************************************************
