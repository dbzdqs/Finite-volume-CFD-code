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
Subroutine GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,E,N1,N2)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,E,P
Intent(Out)::N1,N2

Integer::Dim,P,N1,N2,E,index,I,I1,I2
Integer,Dimension(1:Dim,1:4)::Corn,Neib
!===========================================================================================
if(Corn(E,4) /= 0) then !---------------- NE is a Quadrilateral -----------------
!Part 1:
	do I=1,4
		if(Corn(E,I) == P) then

			index = I

			if(I == 1) then
				I1 = 4
				I2 = I + 1
			elseif(I == 4) then
			    I1 = I - 1
				I2 = 1
			else
				I1 = I - 1
				I2 = I + 1
			endif

			exit

		endif

	end do

	if(I1 < index) then
		if(index - I1 == 1) then
			N1 = Neib(E,I1)
		else
			N1 = Neib(E,index)
		endif
	else
		if(I1 - index == 1) then
			N1 = Neib(E,index)
		else
			N1 = Neib(E,I1)
		endif
	endif

	if(I2 < index) then
		if(index - I2 == 1) then
			N2 = Neib(E,I2)
		else
			N2 = Neib(E,index)
		endif
	else
		if(I2 - index == 1) then
			N2 = Neib(E,index)
		else
			N2 = Neib(E,I2)
		endif
	endif

else !---------------- E is a Triangle -----------------
!Part 2:
	do I=1,3
		if(Corn(E,I) == P) then
			index = I
			exit
		endif
	end do

	if(index == 1) then
		N1 = Neib(E,2)
		N2 = Neib(E,3)
	elseif(index == 2) then
		N1 = Neib(E,3)
		N2 = Neib(E,1)
	elseif(index == 3) then
		N1 = Neib(E,1)
		N2 = Neib(E,2)
	endif

endif
!===========================================================================================
End Subroutine GetTwoNeibourSharingCorner
!*********************************************************************************************
