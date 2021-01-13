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
Subroutine getLargestAngleOfQuad(Dim,Corn,Neib,X,Y,Elm,Ang,P)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,X,Y,Elm
Intent(InOut)::Ang,P

Integer::Dim,Elm,I,J,A,B,C,D,E,index,P
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent,isNotChevron,IsDiagonalInsideQuad,QuadIsChevron
Real(8)::GetAngle,PI,Ang
Real(8),Dimension(1:4)::Angle
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================  
!Part 1:
QuadIsChevron = .Not. isNotChevron(Dim,Corn,X,Y,Elm) 

A = Corn(Elm,1)
B = Corn(Elm,2)
C = Corn(Elm,3)
D = Corn(Elm,4)

if(QuadIsChevron) then
!Part 2:
	if(IsDiagonalInsideQuad(Dim,Corn,X,Y,Elm,A,C) .And. .Not. IsDiagonalInsideQuad(Dim,Corn,X,Y,Elm,B,D)) then
	
		Angle(1) = GetAngle(Dim,A,B,C,X,Y) + GetAngle(Dim,A,C,D,X,Y) !----- Angle at A ----- 
		Angle(2) = GetAngle(Dim,B,A,C,X,Y)							 !----- Angle at B -----
		Angle(3) = GetAngle(Dim,C,A,B,X,Y) + GetAngle(Dim,C,A,D,X,Y) !----- Angle at C -----
		Angle(4) = GetAngle(Dim,D,A,C,X,Y)							 !----- Angle at D -----

	elseif(IsDiagonalInsideQuad(Dim,Corn,X,Y,Elm,B,D) .And. .Not. IsDiagonalInsideQuad(Dim,Corn,X,Y,Elm,A,C)) then
	
		Angle(1) = GetAngle(Dim,A,B,D,X,Y)                           !----- Angle at A ----- 
		Angle(2) = GetAngle(Dim,B,A,D,X,Y) + GetAngle(Dim,B,C,D,X,Y) !----- Angle at B -----
		Angle(3) = GetAngle(Dim,C,B,D,X,Y)                           !----- Angle at C -----
		Angle(4) = GetAngle(Dim,D,A,B,X,Y) + GetAngle(Dim,D,B,C,X,Y) !----- Angle at D -----

	else
	
		Angle(1) = GetAngle(Dim,A,B,D,X,Y) !----- Angle at A -----
		Angle(2) = GetAngle(Dim,B,A,C,X,Y) !----- Angle at B -----
		Angle(3) = GetAngle(Dim,C,B,D,X,Y) !----- Angle at C -----
		Angle(4) = GetAngle(Dim,D,A,C,X,Y) !----- Angle at D -----

	endif
	
else
!Part 3:
	Angle(1) = GetAngle(Dim,A,B,D,X,Y) !----- Angle at A -----
	Angle(2) = GetAngle(Dim,B,A,C,X,Y) !----- Angle at B -----
	Angle(3) = GetAngle(Dim,C,B,D,X,Y) !----- Angle at C -----
	Angle(4) = GetAngle(Dim,D,A,C,X,Y) !----- Angle at D -----	

endif
!Part 4:
Ang = Angle(1)
P = Corn(Elm,1)

do I=2,4

	if(Ang < Angle(I)) then

		Ang = Angle(I)
		P = Corn(Elm,I)

	endif

end do

!===========================================================================================
End Subroutine getLargestAngleOfQuad 
!*********************************************************************************************
