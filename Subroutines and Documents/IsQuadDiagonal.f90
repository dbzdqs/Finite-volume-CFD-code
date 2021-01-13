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
Function IsQuadDiagonal(Dim,Corn,P1,P2,Q)
Implicit None
!===========================================================================================
Intent(In)::P1,P2,Q

Integer::Dim,P1,P2,Q,I,I1,I2
Integer,Dimension(1:Dim,1:4)::Corn
Logical::IsQuadDiagonal,hasCorner
!===========================================================================================
if(hasCorner(Dim,Corn,Q,P1) .And. hasCorner(Dim,Corn,Q,P2)) then
!Part 1:
	do I=1,4
		if(Corn(Q,I) == P1) I1 = I
		if(Corn(Q,I) == P2) I2 = I
	end do
!Part 2:
	if(I2 > I1) then
		if(I2 - I1 == 2) then
			IsQuadDiagonal = .True.
		else
			IsQuadDiagonal = .False.
		endif
	else
		if(I1 - I2 == 2) then
			IsQuadDiagonal = .True.
		else
			IsQuadDiagonal = .False.
		endif
	endif

else
	IsQuadDiagonal = .False.
endif
!===========================================================================================
End Function IsQuadDiagonal 
!*********************************************************************************************
