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
Subroutine FixInvertedElements(Dim,Corn,X,Y,TElms,TEC,QElms,QEC,TEF,QEF,Point,x_old,y_old,DL,Lq_x,Lq_y,Lq)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,QElms,QEC,TElms,TEC,TEF,QEF,Point,x_old,y_old,DL,Lq_x,Lq_y,Lq
Intent(Out)::X,Y

Integer::Dim,QEC,TEC,Point,A,B,C,I
Integer,Dimension(1:1000)::TElms,QElms,TEF,QEF
Integer,Dimension(1:Dim,1:4)::Corn
Logical::ElementInverted,Divergence
Real(8)::DL,L,Lq_x,Lq_y,Lq,x_old,y_old
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
A = 1
B = 2
C = 2
L = A*DL/B

X(Point) = x_old + (Lq_x/Lq)*L
Y(Point) = y_old + (Lq_y/Lq)*L
!Part 2:
if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
	Divergence = .True.
else
	Divergence = .False.
endif

if(Divergence) then
!Part 3:
	do I=1,10
		C = C + 1
		L = Real(DL)/C
		X(Point) = x_old + (Lq_x/Lq)*L
		Y(Point) = y_old + (Lq_y/Lq)*L
        Print *,'1/C:',C,1/C

		if(.Not. ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
			print *,'Exit in:',I
			exit	
		endif

	end do

	if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
		print *,'No Smooth!!!'
		X(Point) = x_old
		Y(Point) = y_old
	endif
else
!Part 4:
	do I=1,10
		A = A + 1
		B = B + 1
		L = Real(A*DL)/B
		X(Point) = x_old + (Lq_x/Lq)*L
		Y(Point) = y_old + (Lq_y/Lq)*L
Print *,'A/B:',Real(A)/B

		if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
			print *,'Exit in:',I
			exit	
		endif

	end do
endif

if(.Not. Divergence) then
print *,'Not Divergence'
	A = A - 1
	B = B - 1
	L = A*DL/B
	X(Point) = x_old + (Lq_x/Lq)*L
	Y(Point) = y_old + (Lq_y/Lq)*L

endif
!===========================================================================================
End Subroutine FixInvertedElements
!*********************************************************************************************
