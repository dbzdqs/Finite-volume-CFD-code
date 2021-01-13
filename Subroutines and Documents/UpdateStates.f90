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
Subroutine UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles) ! Tested
Implicit None
!===========================================================================================
Intent(In)::Dim,X,Y,FrontEdges,Fronts,Corn,Neib
Intent(Out)::States,Angles

Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,I,J,K,index,Fronts,LTC,RTC,corner1,corner2,LV,RV,ME,Prev,C,Next
Real(8)::U1,V1,U2,V2,InnerProduct,Temp,pi,Threshold,Norm1,Norm2
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Integer,Dimension(1:1000)::LeftTris,RightTris
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:Dim,1:2)::Angles
!===========================================================================================
do I=1,Fronts

	if(States(I) /= -1) then
!Part 1:
		Angles(I,LeftVertex) = 0
		Angles(I,RightVertex) = 0

		LV = FrontEdges(I,LeftVertex)
		RV = FrontEdges(I,RightVertex)
		ME = FrontEdges(I,Element)

		!---------------------- Determining triangles sharing common corner on both sides of front edge -----------------------
!Part 2:
		LTC = 1                         ! LTC (Left Triangles Count)
		RTC = 1		                    ! RTC (Right Triangles Count)

		LeftTris(LTC) = ME
		Prev = ME
		C = RV

		do
			do J=1,3
				if(Corn(Prev,J) == C) then
					Next = Neib(Prev,J)
					exit
				endif
			end do
			
			if(Next == 0) then
                exit    
            elseif(Corn(Next,4) /= 0) then
                exit    
            endif
			
			do J=1,3
				if(Corn(Prev,J)/=LV .And. Corn(Prev,J)/=C) then
					C = Corn(Prev,J)
					exit	
				endif
			end do

			LTC = LTC + 1
			LeftTris(LTC) = Next

			Prev = Next
					
		end do
!Part 3:
		RightTris(RTC) = ME
		Prev = ME
		C = LV

		do
			do J=1,3
				if(Corn(Prev,J) == C) then
					Next = Neib(Prev,J)
					exit
				endif
			end do
			
			if(Next == 0) then
                exit    
            elseif(Corn(Next,4) /= 0) then
                exit    
            endif
			
			do J=1,3
				if(Corn(Prev,J)/=RV .And. Corn(Prev,J)/=C) then
					C = Corn(Prev,J)
					exit	
				endif
			end do

			RTC = RTC + 1
			RightTris(RTC) = Next

			Prev = Next
					
		end do
		
		!------------------------------ Calculating Angle on the left side of front edge --------------------------------------
!Part 4:
		do J=1,LTC

			do K=1,3
				if(Corn(LeftTris(J),K) == FrontEdges(I,LeftVertex)) then
					index = K
				endif
			end do

			if(index == 1) then
				corner1 = Corn(LeftTris(J),2)
				corner2 = Corn(LeftTris(J),3)
			elseif(index == 2) then
				corner1 = Corn(LeftTris(J),1)
				corner2 = Corn(LeftTris(J),3)
			elseif(index == 3) then
				corner1 = Corn(LeftTris(J),1)
				corner2 = Corn(LeftTris(J),2)
			endif

			U1=X(corner1) - X(FrontEdges(I,LeftVertex))
			V1=Y(corner1) - Y(FrontEdges(I,LeftVertex))

			U2=X(corner2) - X(FrontEdges(I,LeftVertex))	     
			V2=Y(corner2) - Y(FrontEdges(I,LeftVertex))

			InnerProduct= U1*U2 + V1*V2                       

			Temp=U1*U1 + V1*V1
			Norm1=DSQRT(Temp)
			
			Temp=U2*U2 + V2*V2
			Norm2=DSQRT(Temp)
																   
			Temp=InnerProduct/(Norm1*Norm2)

			if(Temp > 1) then
				Temp = 1.0
			elseif(Temp < -1) then
				Temp = -1.0
			endif

			Angles(I,LeftVertex) = Angles(I,LeftVertex) + ACOS(Temp)

		end do
		!------------------------------ Calculating Angle on the right side of front edge -------------------------------------
!Part 5:
		do J=1,RTC

			do K=1,3
				if(Corn(RightTris(J),K) == FrontEdges(I,RightVertex)) then
					index = K
				endif
			end do

			if(index == 1) then
				corner1 = Corn(RightTris(J),2)
				corner2 = Corn(RightTris(J),3)
			elseif(index == 2) then
				corner1 = Corn(RightTris(J),1)
				corner2 = Corn(RightTris(J),3)
			elseif(index == 3) then
				corner1 = Corn(RightTris(J),1)
				corner2 = Corn(RightTris(J),2)
			endif

			U1=X(corner1) - X(FrontEdges(I,RightVertex))
			V1=Y(corner1) - Y(FrontEdges(I,RightVertex))

			U2=X(corner2) - X(FrontEdges(I,RightVertex))	     
			V2=Y(corner2) - Y(FrontEdges(I,RightVertex))

			InnerProduct= U1*U2 + V1*V2                       

			Temp=U1*U1 + V1*V1
			Norm1=DSQRT(Temp)
			
			Temp=U2*U2 + V2*V2
			Norm2=DSQRT(Temp)
																   
			Temp=InnerProduct/(Norm1*Norm2)
			
			if(Temp > 1) then
				Temp = 1.0
			elseif(Temp < -1) then
				Temp = -1.0
			endif

			Angles(I,RightVertex) = Angles(I,RightVertex) + ACOS(Temp)

		end do		                       
!Part 6:									   
		pi = 3.14159265358979d0
		Threshold = 3*pi/4
									   
		!---------------------------------- Following conditions determines State of Main Front -------------------------------

		if(Angles(I,LeftVertex)>=Threshold .And. Angles(I,RightVertex)>=Threshold) then
			States(I)=0
		elseif(Angles(I,LeftVertex)>=Threshold .And. Angles(I,RightVertex)<Threshold) then
			States(I)=1
		elseif(Angles(I,LeftVertex)<Threshold .And. Angles(I,RightVertex)>=Threshold)	then
			States(I)=2
		elseif(Angles(I,LeftVertex)<Threshold .And. Angles(I,RightVertex)<Threshold) then
			States(I)=3	
		endif

	endif

end do
!===========================================================================================
End Subroutine UpdateStates
!*********************************************************************************************
