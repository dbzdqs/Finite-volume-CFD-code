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
Function SegmentJointIsPossible(Dim,X,Y,A,B,Points,PC)
Implicit None
!===========================================================================================
Intent(In)::Dim,X,Y,A,B,Points,PC

Integer::Dim,PC,A,B,C,D,I
Integer,Dimension(1:1000)::Points
Logical::SegmentJointIsPossible,IntersectionOccur
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
SegmentJointIsPossible = .True.

do I=1,PC
!Part 1:
	C = Points(I)

	if(I < PC) then
	
		D = Points(I+1)

	else
		
		D = Points(1)

	endif
!Part 2:
	if(IntersectionOccur(Dim,A,B,C,D,X,Y)) then
		
		SegmentJointIsPossible = .False.
		exit

	endif 

end do
!===========================================================================================
End Function SegmentJointIsPossible
!*********************************************************************************************
