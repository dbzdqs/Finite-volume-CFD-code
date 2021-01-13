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
Subroutine Smooth(Dim,NC,NBE,BEP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Elm,X,Y)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,NBE,BEP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Elm
Intent(Inout)::X,Y

Integer,Parameter::Processed=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2

Integer::Dim,NC,NBE,Fronts,QEC,TEC,CC,EC,I,J,K,M,LV,RV,Q_Count,Ni,Ni1,Ni2,Nj,Nj1,Nj2,E,Elm,Ii,Ij,Ik,Il,Q1,Q2
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::QElms,TElms,Corners,Elms,TEF,QEF
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::BEP
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::IsOnExteriorBoundary,IsInTheList,ElementInverted,QCorner,hasCorner,isOnTheBoundary,TwoQuadsAreAdjacent,QuadsAreAdjacent
Real(8)::x_old,y_old,GetNorm,DL,pi,LE,SE,L,Tr,Vx,Vy,DeltaAx,DeltaAy,Ld,La,Lq,DeltaBx,DeltaBy,Pxi,Pyi,Pxi1,Pyi1
Real(8)::Pxi2,Pyi2,norm1,norm2,PBx1,PBy1,PBx2,PBy2,Qx,Qy,temp,s,t,DeltaCx,DeltaCy,DeltaIx,DeltaIy,Deltax,Deltay,Dx0,Dy0
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
pi = 3.14159265358979d0

!Part 1:

!---------------------------------- Calculating Tr -----------------------------------------
LV = FrontEdges(1,LeftVertex)
RV = FrontEdges(1,RightVertex)

SE = GetNorm(Dim,LV,RV,X,Y)
LE = SE
	!------------ Finding LargestEdge(LE) and SmallestEdge(SE) on the Boundary -------------
do I=1,Fronts
	if(States(I)/=Processed) then 
		L = GetNorm(Dim,FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex),X,Y)
		if(LE < L) LE = L
		if(SE > L) SE = L
	endif
end do

Tr = LE/SE

Print*,'>>>>>>>>>>>>>>>>Tr:',Tr

! ------------------------------------------- Smoothing ------------------------------------

do I=1,4

	DL = 0 !--------------------- Initialization of DL(Desired Length) ---------------------
	CC = 0
	EC = 0
	Q_Count = 0
	Ni = newQuad(I)
	x_old = X(Ni)
	y_old = Y(Ni)
	Vx = 0
	Vy = 0
    TwoQuadsAreAdjacent = .False.

	if(.Not. IsOnExteriorBoundary(Dim,NBE,BEP,Ni)) then
!Part 2:
		print*,'>>>>> Smoothing Process: <<<<<'
		
        Call GetSurroundingElements(Dim,Corn,Neib,Ni,Elm,TElms,QElms,TEC,QEC)
		Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms,QEC,TElms,TEC,TEF,QEF)
!Part 3:			
		do J=1,QEC
			
			E = QElms(J)
				
			do K=1,4 !--------------- Finding Corresponding Index of Ni ----------------	
				if(Corn(E,K) == Ni) then 
					Ii = K
					exit
				endif	
			end do

			if(Ii == 1) then
				Ij = 2 
				Ik = 3
				Il = 4
			elseif(Ii == 2) then
				Ij = 3 
				Ik = 4
				Il = 1
			elseif(Ii == 3) then
				Ij = 4 
				Ik = 1
				Il = 2
			elseif(Ii == 4) then
				Ij = 1 
				Ik = 2
				Il = 3
			endif

			Vx = Vx + X(Corn(E,Ij)) - X(Corn(E,Ik)) + X(Corn(E,Il))
			Vy = Vy + Y(Corn(E,Ij)) - Y(Corn(E,Ik)) + Y(Corn(E,Il))

		end do
		
		!------------------------------------- Calculating V' -----------------------------------

		Vx = Vx/QEC
		Vy = Vy/QEC

		!------------------------ Amount of change of location of Ni ----------------------------

		DeltaAx = Vx - X(Ni)
		DeltaAy = Vy - Y(Ni)

        if(QEC == 2) then
            Q1 = QElms(1)
			Q2 = QElms(2)
            TwoQuadsAreAdjacent = QuadsAreAdjacent(Dim,Neib,Q1,Q2)
        endif
        
		if(TwoQuadsAreAdjacent) then
!Part 4:
		!--------------------------------- Length Adjustment ------------------------------------

            Q1 = QElms(1)
			Q2 = QElms(2)
            
			do J=1,4 !------------------------------- Finding Nj ------------------------------- 
				do K=1,4
					if(Corn(Q1,J) == Corn(Q2,K) .And. Corn(Q1,J) /= Ni) then
						Nj = Corn(Q1,J)
					endif
				end do
			end do

			!-------------------- Finding Ni-1 -----------------------
			do J=1,4
				if(Corn(Q1,J) == Ni) then
					Ii = J
					exit
				endif
			end do

			if(Ii==1) then
				if(Corn(Q1,2)/=Nj) then
					Ni1 = Corn(Q1,2)	
				else
					Ni1 = Corn(Q1,4)
				endif
			elseif(Ii==4) then
				if(Corn(Q1,1)/=Nj) then
					Ni1 = Corn(Q1,1)	
				else
					Ni1 = Corn(Q1,3)
				endif
			else
				if(Corn(Q1,Ii+1)/=Nj) then
					Ni1 = Corn(Q1,Ii+1)	
				else
					Ni1 = Corn(Q1,Ii-1)
				endif
			endif 

			!-------------------- Finding Ni+1 -----------------------
			do J=1,4
				if(Corn(Q2,J) == Ni) then
					Ii = J
					exit
				endif
			end do

			if(Ii==1) then
				if(Corn(Q2,2)/=Nj) then
					Ni2 = Corn(Q2,2)	
				else
					Ni2 = Corn(Q2,4)
				endif
			elseif(Ii==4) then
				if(Corn(Q2,1)/=Nj) then
					Ni2 = Corn(Q2,1)	
				else
					Ni2 = Corn(Q2,3)
				endif
			else
				if(Corn(Q2,Ii+1)/=Nj) then
					Ni2 = Corn(Q2,Ii+1)	
				else
					Ni2 = Corn(Q2,Ii-1)
				endif
			endif

			!--------------- Finding Connecting Corners to Nj (Nj-1 and Nj+1)-------------------
			do J=1,4 
				if(Corn(Q1,J)/=Ni .And. Corn(Q1,J)/=Ni1 .And. Corn(Q1,J)/=Nj) then
					Nj1 = Corn(Q1,J)
					exit
				endif
			end do

			do J=1,4 
				if(Corn(Q2,J)/=Ni .And. Corn(Q2,J)/=Ni2 .And. Corn(Q2,J)/=Nj) then
					Nj2 = Corn(Q2,J)
					exit
				endif
			end do

			!--------------------------------- Calculating Ld ----------------------------------

            Ld = (GetNorm(Dim,Nj,Nj1,X,Y) + GetNorm(Dim,Nj,Nj2,X,Y))/2 
			
			!--------------------------------- Calculating La ----------------------------------

			La = DSQRT(((Vx-X(Nj))*(Vx-X(Nj))) + ((Vy-Y(Nj))*(Vy-Y(Nj))))

			!---------------------- Amount of change in position of Ni -------------------------

			DeltaBx = X(Nj) - X(Ni) + (DeltaAx + X(Ni) - X(Nj))*(Ld/La)
			DeltaBy = Y(Nj) - Y(Ni) + (DeltaAy + Y(Ni) - Y(Nj))*(Ld/La)

		!-------------------------------  Angular Smoothness -----------------------------------

			Pxi = X(Ni) - X(Nj) !------------------ Calculating Vector Pi ----------------------
			Pyi = Y(Ni) - Y(Nj)

			Pxi1 = X(Ni1) - X(Nj) !---------------- Calculating Vector Pi-1 --------------------
			Pyi1 = Y(Ni1) - Y(Nj)

			Pxi2 = X(Ni2) - X(Nj) !---------------- Calculating Vector Pi+1 --------------------
			Pyi2 = Y(Ni2) - Y(Nj)

			!----------------------- Calculating Bisector Vector PB1 ---------------------------

			norm1 = DSQRT(Pxi1*Pxi1 + Pyi1*Pyi1) !-------------- Norm of Pi-1 ------------------
			norm2 = DSQRT(Pxi2*Pxi2 + Pyi2*Pyi2) !-------------- Norm of Pi+1 ------------------

			PBx1 = norm1*Pxi2 + norm2*Pxi1
			PBy1 = norm1*Pyi2 + norm2*Pyi1

			!----------------------- Calculating Bisector Vector PB2 ---------------------------

			norm1 = DSQRT(Pxi*Pxi + Pyi*Pyi) !------------------ Norm of Pi --------------------
			norm2 = DSQRT(PBx1*PBx1 + PBy1*PBy1)!--------------- Norm of PB1 -------------------
			
			PBx2 = norm1*PBx1 + norm2*Pxi
			PBy2 = norm1*PBy1 + norm2*Pyi

			!--------- Calculating Coordinates of point Q based on Eberly(2000) Method ---------

			Dx0 = X(Ni2) - X(Ni1)
			Dy0 = Y(Ni2) - Y(Ni1)

			!----------------------------->>> Notice: Here D1 = PB2 <<<-------------------------

			Deltax = X(Nj) - X(Ni1)
			Deltay = Y(Nj) - Y(Ni1)

			temp = Dx0*PBy2 - PBx2*Dy0

			if(temp /= 0) then
				
				s = (Deltax*PBy2 - PBx2*Deltay)/temp
				t =	(Deltax*Dy0 - Dx0*Deltay)/temp

			endif

			Qx = X(Ni1) + s*Dx0
			Qy = Y(Ni1) + s*Dy0

			!--------------------------------- Calculating Lq ----------------------------------

			Lq = DSQRT((Qx - X(Nj))*(Qx - X(Nj)) + (Qy - Y(Nj))*(Qy - Y(Nj)))

			!----------------------------- Determining Length of PB2 ---------------------------

			if(Ld > Lq) then
				temp = (Lq + Ld)/2
			else
				temp = Ld
			endif
			
			norm1 = DSQRT(PBx2*PBx2 + PBy2*PBy2)
			
			PBx2 = PBx2*temp/norm1
			PBy2 = PBy2*temp/norm1
			
			DeltaCx = PBx2 - Pxi
			DeltaCy = PBy2 - Pyi
			
			DeltaIx = (DeltaBx + DeltaCx)/2
			DeltaIy = (DeltaBy + DeltaCy)/2 

		endif

		if(Tr < 2.5) then !---- No Modification in Smoothing Method Proposed by Blacker (1991) -----

!Part 5:
			
			if(TwoQuadsAreAdjacent) then
				X(Ni) = X(Ni) + DeltaIx
				Y(Ni) = Y(Ni) + DeltaIy
			else
				X(Ni) = X(Ni) + DeltaAx
				Y(Ni) = Y(Ni) + DeltaAy
			endif

			!--------------------- Finding non-Quad Corners surrounding Ni -------------------------
			do J=1,TEC
				E = TElms(J)
				do K=1,3
					QCorner = .False.
					do M=1,QEC
						if(hasCorner(Dim,Corn,QElms(M),Corn(E,K))) then
							QCorner = .True.
							exit
						endif	
					end do

					if(Corn(E,K)/=Ni .And. (.Not. QCorner) .And. (.Not. IsInTheList(Corners,CC,Corn(E,K)))) then
						CC = CC + 1
						Corners(CC) = Corn(E,K)
						EC = EC + 1
						Elms(EC) = E
					endif
				end do
			end do

        elseif(Tr >= 2.5 .And. Tr <= 20) then !-- Using Ld Proposed by Owen 1999 (Main Reference) --

!Part 6:
			
			!--------------------- Finding non-Quad Corners surrounding Ni -------------------------
			do J=1,TEC
				E = TElms(J)
				do K=1,3
					QCorner = .False.
					do M=1,QEC
						if(hasCorner(Dim,Corn,QElms(M),Corn(E,K))) then
							QCorner = .True.
							exit
						endif	
					end do

					if(Corn(E,K)/=Ni .And. (.Not. QCorner) .And. (.Not. IsInTheList(Corners,CC,Corn(E,K)))) then
						CC = CC + 1
						Corners(CC) = Corn(E,K)
						EC = EC + 1
						Elms(EC) = E
					endif
				end do
			end do

			if(TwoQuadsAreAdjacent) then
				
				Ld = 0

				do J=1,CC
					Ld = Ld + GetNorm(Dim,Ni,Corners(J),X,Y) 	
				end do
				
				Ld = Ld + GetNorm(Dim,Ni1,Nj1,X,Y)
				Ld = Ld + GetNorm(Dim,Ni2,Nj2,X,Y)
				Ld = Ld + GetNorm(Dim,Nj,Nj1,X,Y)
				Ld = Ld + GetNorm(Dim,Nj,Nj2,X,Y)
				
				Ld = Ld/(4 + CC) 

				X(Ni) = X(Ni) + DeltaIx
				Y(Ni) = Y(Ni) + DeltaIy

				Pxi = X(Ni)-X(Nj)
				Pyi = Y(Ni)-Y(Nj) 
				norm1 = DSQRT(Pxi*Pxi + Pyi*Pyi)

				X(Ni) = X(Nj) + Pxi*Ld/norm1
				Y(Ni) = Y(Nj) + Pyi*Ld/norm1
				
			else
				X(Ni) = X(Ni) + DeltaAx
				Y(Ni) = Y(Ni) + DeltaAy
			endif

		elseif(Tr > 20) then

!Part 7:
			!--------------------- Finding non-Quad Corners surrounding Ni -------------------------
			do J=1,TEC
				E = TElms(J)
				do K=1,3
					QCorner = .False.
					do M=1,QEC
						if(hasCorner(Dim,Corn,QElms(M),Corn(E,K))) then
							QCorner = .True.
							exit
						endif	
					end do

					if(Corn(E,K)/=Ni .And. (.Not. QCorner) .And. (.Not. IsInTheList(Corners,CC,Corn(E,K)))) then
						CC = CC + 1
						Corners(CC) = Corn(E,K)
						EC = EC + 1
						Elms(EC) = E
					endif
				end do
			end do

			if(TwoQuadsAreAdjacent) then

				Ld = 0

				do J=1,CC
					Ld = Ld + GetNorm(Dim,Ni,Corners(J),X,Y) 
				end do

				Ld = Ld + GetNorm(Dim,Ni,Ni1,X,Y)
				Ld = Ld + GetNorm(Dim,Ni,Ni2,X,Y)
				Ld = Ld + GetNorm(Dim,Ni,Nj,X,Y)
				Ld = Ld/(3+CC)

				X(Ni) = X(Ni) + DeltaIx
				Y(Ni) = Y(Ni) + DeltaIy

				Pxi = X(Ni)-X(Nj)
				Pyi = Y(Ni)-Y(Nj) 
				norm1 = DSQRT(Pxi*Pxi + Pyi*Pyi)

				X(Ni) = X(Nj) + Pxi*Ld/norm1
				Y(Ni) = Y(Nj) + Pyi*Ld/norm1

			else

				X(Ni) = X(Ni) + DeltaAx
				Y(Ni) = Y(Ni) + DeltaAy

			endif

		endif

		!----------------------------- Fix Inverted Elements if Exist ------------------------------

		Pxi = X(Ni) - x_old
		Pyi = Y(Ni) - y_old
		Ld = DSQRT(Pxi*Pxi + Pyi*Pyi)
!Part 8:
		if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then

			Call FixInvertedElements(Dim,Corn,X,Y,TElms,TEC,QElms,QEC,TEF,QEF,Ni,x_old,y_old,Ld,Pxi,Pyi,Ld)
            
		endif
!Part 9:
		!------------------- Laplacian Smooth of non-Quad Corners surrounding Ni ---------------
		do J=1,CC
			if(.Not. isOnTheBoundary(Dim,Fronts,Corners(J),FrontEdges,States)) then
				print *,'Laplacian Smooth:',Elms(J),Corners(J),':',Corn(Elms(J),1),Corn(Elms(J),2),Corn(Elms(J),3)
				Call LaplacianSmooth(Dim,Corners(J),Elms(J),Corn,Neib,X,Y)
			endif
		end do

	endif

end do
	
!===========================================================================================
End Subroutine Smooth
!*********************************************************************************************
