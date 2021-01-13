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
Subroutine GetNextTriangle(Dim,NC,Corn,Ti,Ei,Ti1)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,Corn,Ti,Ei
Intent(Out)::Ti1

Integer::Dim,NC,Ti,Ti1,I,J
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:Dim,1:4)::Corn
!===========================================================================================
!Part 1:
do I=1,NC

	do J=1,3
		if(Corn(I,J) == Ei(1) .And. I /= Ti .And. Corn(I,4) == 0) then
			if(J == 1) then
				if(Corn(I,2) == Ei(2) .Or. Corn(I,3) == Ei(2)) then
					Ti1 = I	
				endif
			elseif(J == 2) then
				if(Corn(I,1) == Ei(2) .Or. Corn(I,3) == Ei(2)) then
					Ti1 = I	
				endif
			elseif(J == 3) then
				if(Corn(I,1) == Ei(2) .Or. Corn(I,2) == Ei(2)) then
					Ti1 = I	
				endif
			endif
		endif
	end do

end do
!===========================================================================================
End Subroutine GetNextTriangle
!*********************************************************************************************
