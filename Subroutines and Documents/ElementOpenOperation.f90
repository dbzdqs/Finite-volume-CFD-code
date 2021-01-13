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
Subroutine ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,C,A,B,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,ME,NBE,BFP,C,A,B
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,ME,NE,I,J,K,TEC,QEC,Elm,A,Ai,B,C,Ci,Na,P,getNeibour,PC1,PC2,EC1,EC2,ABElm1,ABElm2,newElement,A2
Integer,Dimension(1:1000)::TElms,QElmsA,TEF,QEF,PList1,PList2,EList1,EList2
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,areAdjacent,ElementInverted,isQuadNode,QuadIsInverted
Real(8)::x_old_A,y_old_A,x_A1,y_A1,x_A2,y_A2
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
done = .False.
!Part 1:
x_old_A = X(A)
y_old_A = Y(A)

Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElmsA,TEC,Na)
Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElmsA,Na,TElms,TEC,TEF,QEF)

PC1 = 0 !------------ Counter of First List of Points (PList1) -----------------
PC2 = 0	!------------ Counter of Second List of Points (PList2) ----------------

EC1 = 0 !------------ Counter of First List of Elements (EList1) ---------------
EC2 = 0	!------------ Counter of Second List of Elements (EList2) --------------

!Part 2:

!--------------------- Finding Indices of A and C in ME ------------------------

do I=1,4

	if(Corn(ME,I) == A) Ai = I
	if(Corn(ME,I) == C) Ci = I

end do

!---------------- Finding NE (Element Neibour of ME on AC side) ----------------
				
NE = getNeibour(Dim,Corn,Neib,A,C,ME)

!--------------------- Adding 'C' corner to both Point Lists -------------------

Call addToList(PList1,PC1,C)
Call addToList(PList2,PC2,C)

!Part 3:

!------------------------ Processing PList1 and EList1 -------------------------

P = C
Elm = ME

do

	do J=1,4

		if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P) then

			P = Corn(Elm,J)

			Call addToList(PList1,PC1,P)
			Call addToList(EList1,EC1,Elm)

			if(P /= B) then

				Elm = getNeibour(Dim,Corn,Neib,A,P,Elm)

			endif

			exit

		endif
	
	end do 
	
	if(P == B) then

		ABElm1 = Elm !---------- Element in Elist1 having AB Edge ----------
		exit

	endif	
  
end do !--------------- End of Processing PList1 and EList1 ----------------

!---------------------- Processing PList2 and EList2 -----------------------

P = C
Elm = NE

do

	do J=1,4

		if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P) then

			P = Corn(Elm,J)

			Call addToList(PList2,PC2,P)
			Call addToList(EList2,EC2,Elm)

			if(P /= B) then

				Elm = getNeibour(Dim,Corn,Neib,A,P,Elm)
			
			endif

			exit

		endif
	
	end do 
	
	if(P == B) then
	
		ABElm2 = Elm !---------- Element in Elist2 having AB Edge ----------
		exit

	endif 	
  
end do !--------------- End of Processing PList2 and EList2 ----------------

!Part 4:
	 
Call CenterOfPolygon(Dim,PC1,PList1,X,Y,x_A1,y_A1)
Call CenterOfPolygon(Dim,PC2,PList2,X,Y,x_A2,y_A2)

!Part 5:

!--->>> Notice: Only one point is added to mesh and A is mapped to A1 coordinates

X(A) = x_A1
Y(A) = y_A1 

NP = NP + 1
A2 = NP

X(A2) = x_A2
Y(A2) = y_A2

do J=1,EC2
	
	Elm = EList2(J)

	do K=1,4

		if(Corn(Elm,K) == A) then
			
			Corn(Elm,K) = A2
			exit

		endif

	end do

end do

!Part 6:

NC = NC + 1
newElement = NC

!--->>> Notice: We use ME corner indices to set newElement corners in order to KEEP newElement CCW.

Corn(newElement,Ai) = B
Corn(newElement,Ci) = A2

do J=1,4

	if(areAdjacent(Dim,Corn,A,Corn(ME,J),ME) .And. Corn(ME,J) /= C) then

		Corn(newElement,J) = A 

	endif

	if(areAdjacent(Dim,Corn,C,Corn(ME,J),ME) .And. Corn(ME,J) /= A) then

		Corn(newElement,J) = C 

	endif

end do

!-------------------- Determining newElement Neibours ------------------

Call setNeibour(Dim,Corn,Neib,newElement,A,B,ABElm1)
Call setNeibour(Dim,Corn,Neib,newElement,A2,B,ABElm2)
Call setNeibour(Dim,Corn,Neib,newElement,A,C,ME)
Call setNeibour(Dim,Corn,Neib,newElement,A2,C,NE)

!------------------- Determining Neibours of Neibours ------------------

Call setNeibour(Dim,Corn,Neib,ABElm1,A,B,newElement)
Call setNeibour(Dim,Corn,Neib,ABElm2,A2,B,newElement)
Call setNeibour(Dim,Corn,Neib,ME,A,C,newElement)
Call setNeibour(Dim,Corn,Neib,NE,A2,C,newElement)

!Part 7:

if(ElementInverted(Dim,Corn,X,Y,TElms,QElmsA,TEC,Na,TEF,QEF) .Or. QuadIsInverted(Dim,Corn,X,Y,newElement)) then

	X(A) = x_old_A 
	Y(A) = y_old_A
	
	do J=1,EC2
	
		Elm = EList2(J)

		do K=1,4

			if(Corn(Elm,K) == A2) then
				
				Corn(Elm,K) = A
				exit

			endif

		end do

	end do
	
	!------------------- Restoring Neibours of Neibours ----------------

	Call setNeibour(Dim,Corn,Neib,ABElm1,A,B,ABElm2)
	Call setNeibour(Dim,Corn,Neib,ABElm2,A,B,ABElm1)
	Call setNeibour(Dim,Corn,Neib,ME,A,C,NE)
	Call setNeibour(Dim,Corn,Neib,NE,A,C,ME)
	
	!----------------------- Deleting newElement -----------------------
	
	do J=2,4
	
		Corn(newElement,J) = Corn(newElement,1)

	end do
	
	NC = NC - 1 !------>>>>>> Notice: No Element Added! <<<<<<<<-------- 

else
!Part 8:
	done = .True.

	!------------------------ Smoothing Elements -----------------------

	do J=1,4

		Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(newElement,J),newElement)

	end do

	do J=1,Na
	
		Elm = QElmsA(J) 

		do K=1,4

			Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Elm,K),Elm)

		end do

	end do

endif
 
!===========================================================================================
End Subroutine ElementOpenOperation
!*********************************************************************************************
