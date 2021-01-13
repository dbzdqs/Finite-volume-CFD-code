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
Function hasTwoCoincidentNode(Dim,Corn,X,Y,Elm,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Elm,COINCIDENT_TOLERANCE

Integer::Dim,Elm,A,B,I,index
Integer,Dimension(1:Dim,1:4)::Corn
Logical::hasTwoCoincidentNode
Real(8)::COINCIDENT_TOLERANCE,GetNorm
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

hasTwoCoincidentNode = .False.
!Part 1:
do I=1,4

	A = Corn(Elm,I)

	Call getNextCorner(Dim,Corn,Elm,A,B,index)

	if(GetNorm(Dim,A,B,X,Y) < COINCIDENT_TOLERANCE) then

		hasTwoCoincidentNode = .True.
		exit

	endif

end do
!===========================================================================================
End Function hasTwoCoincidentNode
!*********************************************************************************************
