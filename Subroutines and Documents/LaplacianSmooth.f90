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
Subroutine LaplacianSmooth(Dim,point,E,Corn,Neib,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,point,E,Corn,Neib
Intent(Inout)::X,Y

Integer::Dim,point,E,SC_Count,TEC,QEC,I,J,K,N,index
Integer,Dimension(1:1000)::SC,TElms,QElms,TEF,QEF
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsInTheList,ElementInverted
Real(8)::x_temp,y_temp,x_old,y_old
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
SC_Count = 0 ! Surrounding Corners
TEC = 0
QEC = 0
x_temp = 0
y_temp = 0
x_old = X(point)
y_old = Y(point)
!Part 1:
if(Corn(E,4) == 0) then

	do I=1,3
		if(Corn(E,I) == point) then
			index = I
			exit
		endif
	end do

	if(index-1 < 1) then
		J = 3
	else
		J = index - 1
	endif

	if(index+1 > 3) then
		K = 1
	else
		K = index + 1
	endif

	if(Neib(E,J) /= 0) then
		N = Neib(E,J) 
	elseif(Neib(E,K) /= 0) then
		N = Neib(E,K)
	endif

else

	do I=1,4
		if(Corn(E,I) == point) then
			index = I
			exit
		endif
	end do

	if(index-1 < 1) then
		J = 4
	else
		J = index - 1
	endif

	if(index+1 > 4) then
		K = 1
	else
		K = index + 1
	endif

	if(Neib(E,J) /= 0) then
		N = Neib(E,J) 
	elseif(Neib(E,K) /= 0) then
		N = Neib(E,K)
	endif

endif
!Part 2:
Call GetSurroundingElements(Dim,Corn,Neib,point,N,TElms,QElms,TEC,QEC)
Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms,QEC,TElms,TEC,TEF,QEF)
!Part 3:
do I=1,TEC
	N = TElms(I)
	do J=1,3
		if(Corn(N,J) == point .And. Corn(N,1)/=Corn(N,2)) then
			if(J == 1) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,2))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,2)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,3))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,3)
				endif	
			elseif(J == 2) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,1))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,1)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,3))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,3)
				endif
			elseif(J == 3) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,1))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,1)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,2))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,2)
				endif
			endif
		endif
	end do
end do

do I=1,QEC
	N = QElms(I)
	do J=1,4
		if(Corn(N,J) == point .And. Corn(N,1)/=Corn(N,2)) then
			if(J == 1) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,2))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,2)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,4))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,4)
				endif	
			elseif(J == 2) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,1))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,1)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,3))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,3)
				endif
			elseif(J == 3) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,4))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,4)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,2))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,2)
				endif
			elseif(J == 4) then
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,1))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,1)	
				endif
				if(.Not. IsInTheList(SC,SC_Count,Corn(N,3))) then
					SC_Count = SC_Count + 1
					SC(SC_Count) = Corn(N,3)
				endif
			endif
		endif
	end do
end do
!Part 4:
do I=1,SC_Count
	x_temp = x_temp + X(SC(I))
	y_temp = y_temp + Y(SC(I))
end do

X(point) = x_temp/SC_Count
Y(point) = y_temp/SC_Count
!Part 5:
if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
	X(point) = x_old
	Y(point) = y_old
endif
!===========================================================================================
End Subroutine LaplacianSmooth
!*********************************************************************************************
