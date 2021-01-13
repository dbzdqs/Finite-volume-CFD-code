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
Subroutine ElementOpenningOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,C,A,M,E1,E2,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,C,A,M,E1,E2
Intent(Out)::done
Intent(InOut)::NC,NP,Corn,Neib,X,Y

Integer::Dim,NC,NP,NBE,E1,E2,I,J,K,TEC,QEC,Elm,A,Ai,B,Bi,C,Ci,D,Di,M,P1,getNeibour,EC1,EC2,newElement,A1,A2,N1,N2
Integer,Dimension(1:1000)::TElms,QElms,TEF,QEF,PList1,PList2,EList1,EList2
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent,QuadIsInverted,CCW,AnyElementInverted,IsOnExteriorBoundary,done
Real(8)::DELTA,COINCIDENT_TOLERANCE,temp,SE,GetNorm,ACx,ACy,Vx,Vy,Vnorm
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:

Call CalculateConstants(Dim,NBE,BFP,X,Y,DELTA,COINCIDENT_TOLERANCE)

Call GetSurroundingElements(Dim,Corn,Neib,A,E1,TElms,QElms,TEC,QEC)

!Part 2:

!------------ Finding corners and indices of E1 element ----------------------
do I=1,4
    if(Corn(E1,I) == A) Ai = I
    if(Corn(E1,I) == C) Ci = I
    if(areAdjacent(Dim,Corn,A,Corn(E1,I),E1) .And. Corn(E1,I) /= C) then
        B = Corn(E1,I)
        Bi = I
    endif
    if(areAdjacent(Dim,Corn,C,Corn(E1,I),E1) .And. Corn(E1,I) /= A) then
        D = Corn(E1,I)
        Di = I
    endif
end do

!Part 3:

!---------------- Checking which side E1 and E2 has located -----------------

if(ABS(Ci - Ai) == 1) then

    if(Ci > Ai) then
        CCW = .True.    
    else
        CCW = .False.    
    endif
    
else

    if(Ci > Ai) then
        CCW = .False.    
    else
        CCW = .True.    
    endif
    
endif

N1 = getNeibour(Dim,Corn,Neib,A,C,E1)
N2 = getNeibour(Dim,Corn,Neib,A,M,E2)

!Part 4:

!-------------------- Finding shortest edge connected to 'A' corner -------------------------
SE = GetNorm(Dim,A,C,X,Y)
P1 = C
Elm = E1

do
            
	do J=1,4

		if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P1) then
							
			P1 = Corn(Elm,J)
            Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)
            temp = GetNorm(Dim,A,P1,X,Y) 
            
            if(temp < SE) then
                SE = temp    
            endif
            
			exit

		endif

	end do

	if(P1 == C) exit

end do

SE = SE/2 !------------------- Half size of smallest edge ----------------------------

!Part 5:

!-------------------- Calculating Perpendicular Vector to AC -------------------------

ACx = X(C) - X(A)
ACy = Y(C) - Y(A)

Vx = (-1)*ACy
Vy = ACx

Vnorm = DSQRT(Vx*Vx + Vy*Vy)

Vx = (Vx/Vnorm)*SE !--------- Amount of movement of 'A' corner in both sides --------
Vy = (Vy/Vnorm)*SE

NP = NP + 1
A1 = NP            !-------- new corner A1 is at the direction of prep vector -------

X(A1) = X(A) + Vx
Y(A1) = Y(A) + Vy

NP = NP + 1
A2 = NP            !----- new corner A2 is at the oppo direction of prep vector -----

X(A2) = X(A) - Vx
Y(A2) = Y(A) - Vy

!Part 6:

!--------------------------- Finding EList1 of elements -----------------------------

EC1 = 0 !------------ Counter of First List of Elements (EList1) ---------------
EC2 = 0	!------------ Counter of Second List of Elements (EList2) --------------

P1 = C
Elm = E1

do

	do J=1,4

		if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P1) then

			P1 = Corn(Elm,J)

			Call addToList(EList1,EC1,Elm)

			if(P1 /= M) then

				Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)

			endif

			exit

		endif
	
	end do 
	
	if(P1 == M) then
		exit
	endif	
  
end do !------------------- End of Processing EList1 --------------------

!-------------------------------- Finding EList2 ------------------------

P1 = C
Elm = N1

do

	do J=1,4

		if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P1) then

			P1 = Corn(Elm,J)

			Call addToList(EList2,EC2,Elm)

			if(P1 /= M) then

				Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)
			
			endif

			exit

		endif
	
	end do 
	
	if(P1 == M) then
	
		exit

	endif 	
  
end do !--------------- End of Processing PList2 and EList2 ----------------

!Part 7:

!---------------------------- Adding new element ---------------------------

NC = NC + 1
newElement = NC
Corn(newElement,Di) = C
Corn(newElement,Ai) = M

if(CCW) then
    
    Corn(newElement,Ci) = A2
    Corn(newElement,Bi) = A1
    
    Call setNeibour(Dim,Corn,Neib,newElement,A1,C,E1)
    Call setNeibour(Dim,Corn,Neib,newElement,A1,M,E2)
    Call setNeibour(Dim,Corn,Neib,newElement,A2,C,N1)
    Call setNeibour(Dim,Corn,Neib,newElement,A2,M,N2)
    
    do I=1,EC1
        Elm = EList1(I)
        do J=1,4
            if(Corn(Elm,J) == A) then
                Corn(Elm,J) = A1
                exit
            endif
        end do
    end do
    
    do I=1,EC2
        Elm = EList2(I)
        do J=1,4
            if(Corn(Elm,J) == A) then
                Corn(Elm,J) = A2
                exit
            endif
        end do
    end do
    
    Call setNeibour(Dim,Corn,Neib,E1,A1,C,newElement)
    Call setNeibour(Dim,Corn,Neib,E2,A1,M,newElement)
    Call setNeibour(Dim,Corn,Neib,N1,A2,C,newElement)
    Call setNeibour(Dim,Corn,Neib,N2,A2,M,newElement)
    
else
    
    Corn(newElement,Ci) = A1
    Corn(newElement,Bi) = A2
    
    Call setNeibour(Dim,Corn,Neib,newElement,A1,C,N1)
    Call setNeibour(Dim,Corn,Neib,newElement,A1,M,N2)
    Call setNeibour(Dim,Corn,Neib,newElement,A2,C,E1)
    Call setNeibour(Dim,Corn,Neib,newElement,A2,M,E2)
    
    do I=1,EC1
        Elm = EList1(I)
        do J=1,4
            if(Corn(Elm,J) == A) then
                Corn(Elm,J) = A2
                exit
            endif
        end do
    end do
    
    do I=1,EC2
        Elm = EList2(I)
        do J=1,4
            if(Corn(Elm,J) == A) then
                Corn(Elm,J) = A1
                exit
            endif
        end do
    end do
    
    Call setNeibour(Dim,Corn,Neib,N1,A1,C,newElement)
    Call setNeibour(Dim,Corn,Neib,N2,A1,M,newElement)
    Call setNeibour(Dim,Corn,Neib,E1,A2,C,newElement)
    Call setNeibour(Dim,Corn,Neib,E2,A2,M,newElement)
    
endif

!Part 8:

do I=1,10
    
    if(AnyElementInverted(Dim,Corn,X,Y,TEC,TElms,QEC,QElms)) then !-------------------------- Elements inverted -------------------------------
        
        print*,'Element Inversion after Element Open!',I
        
        SE = SE/2
        Vnorm = DSQRT(Vx*Vx + Vy*Vy)

        Vx = (Vx/Vnorm)*SE !--------- Amount of movement of 'A' corner in both sides --------
        Vy = (Vy/Vnorm)*SE
        
        X(A1) = X(A) + Vx
        Y(A1) = Y(A) + Vy

        X(A2) = X(A) - Vx
        Y(A2) = Y(A) - Vy
        
    else
        
        exit

    endif

end do

if(AnyElementInverted(Dim,Corn,X,Y,TEC,TElms,QEC,QElms)) then

    !---------------------- Undoing open operation ------------------------------
    
    Call DeleteCell(Dim,newElement,Corn)
    
    NC = NC - 1
    
    do I=1,EC1
        Elm = EList1(I)
        do J=1,4
            if(Corn(Elm,J) == A1 .Or. Corn(Elm,J) == A2) then
                Corn(Elm,J) = A
                exit
            endif
        end do
    end do
    
    do I=1,EC2
        Elm = EList2(I)
        do J=1,4
            if(Corn(Elm,J) == A1 .Or. Corn(Elm,J) == A2) then
                Corn(Elm,J) = A
                exit
            endif
        end do
    end do
    
    Call setNeibour(Dim,Corn,Neib,E1,A,C,N1)
    Call setNeibour(Dim,Corn,Neib,N1,A,C,E1)
    Call setNeibour(Dim,Corn,Neib,E2,A,M,N2)
    Call setNeibour(Dim,Corn,Neib,N2,A,M,E2)
    
    done = .False.
    
else

    !----------------------------------------- Smoothing ---------------------------------------
    do I=1,QEC
    
        Elm = QElms(I)
    
        do J=1,4

            if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,Corn(Elm,J))) then
            
		        Call CLS(Dim,Corn,Neib,X,Y,Elm,Corn(Elm,J),COINCIDENT_TOLERANCE)
        
            endif

        end do
    
    end do
    
    done = .True.
    
endif

!===========================================================================================
End Subroutine ElementOpenningOperation
!*********************************************************************************************
