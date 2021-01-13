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
Subroutine NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,E,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,NBE,BFP,E
Intent(InOut)::Corn,Neib,X,Y
Intent(Out)::done

Integer::Dim,NC,NBE,TEC,QEC,E,I,J,C1,C_i1,L,L_i1,R,R_i1,NL,NR,C2,C_i2,L_i2,R_i2,newNL,newNR,N
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsOnExteriorBoundary,done
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

done = .False.

do J=1,4
!Part 1:
	Call GetSurroundingElements(Dim,Corn,Neib,Corn(E,J),E,TElms,QElms,TEC,QEC)

	if(QEC == 2 .And. TEC == 0 .And. .Not. IsOnExteriorBoundary(Dim,NBE,BFP,Corn(E,J))) then

		done = .True.

		C1 = Corn(E,J)
		C_i1 = J
!Part 2:					
		if(C_i1 == 1) then
			L_i1 = 2			!------- L_i(Next Index)
			R_i1 = 4			!------- R_i(Previous Index)
			L = Corn(E,L_i1)
			R = Corn(E,R_i1)	
			NL = Neib(E,1)	    !------- NL (Neibour on the Left side of node 'C')
			NR = Neib(E,4)	    !------- NR (Neibour on the Right side of node 'C')
		elseif(C_i1 == 4) then
			L_i1 = 1
			R_i1 = 3
			L = Corn(E,L_i1)
			R = Corn(E,R_i1)
			NL = Neib(E,4)
			NR = Neib(E,3)	
		else
			L_i1 = C_i1 + 1 
			R_i1 = C_i1 - 1
			L = Corn(E,L_i1)
			R = Corn(E,R_i1)
			NL = Neib(E,C_i1)
			NR = Neib(E,C_i1 - 1) 	
		endif 
!Part 3:					
		do
			if(NL /= NR) exit

			N = NL

			do I=1,4
				if(Corn(N,I)/=C1 .And. Corn(N,I)/=L .And. Corn(N,I)/=R) then
					C2 = Corn(N,I)
					C_i2 = I
				endif
				if(Corn(N,I) == L) L_i2 = I
				if(Corn(N,I) == R) R_i2 = I
			end do

			!------- Specifying new NL ----------
			if(C_i2 > L_i2) then
				if(C_i2 - L_i2 == 1) then
					newNL = Neib(N,L_i2)
				else
					newNL = Neib(N,C_i2)
				endif
			else
				if(L_i2 - C_i2 == 1) then
					newNL = Neib(N,C_i2)
				else
					newNL = Neib(N,L_i2)
				endif
			endif
			!------- Specifying new NR ----------
			if(C_i2 > R_i2) then
				if(C_i2 - R_i2 == 1) then
					newNR = Neib(N,R_i2)
				else
					newNR = Neib(N,C_i2)
				endif
			else
				if(R_i2 - C_i2 == 1) then
					newNR = Neib(N,C_i2)
				else
					newNR = Neib(N,R_i2)
				endif
			endif

			!------------------ Modifying E neibours ------------------

			if(C_i1 > L_i1) then
				if(C_i1 - L_i1 == 1) then
					Neib(E,L_i1) = newNL
				else
					Neib(E,C_i1) = newNL
				endif
			else
				if(L_i1 - C_i1 == 1) then
					Neib(E,C_i1) = newNL
				else
					Neib(E,L_i1) = newNL
				endif
			endif

			if(C_i1 > R_i1) then
				if(C_i1 - R_i1 == 1) then
					Neib(E,R_i1) = newNR
				else
					Neib(E,C_i1) = newNR
				endif
			else
				if(R_i1 - C_i1 == 1) then
					Neib(E,C_i1) = newNR
				else
					Neib(E,R_i1) = newNR
				endif
			endif

			NL = newNL
			NR = newNR

			!-------------- Modify newNL and newNR Neibours ------------
			do I=1,4
				if(Neib(newNL,I) == N) Neib(newNL,I) = E
				if(Neib(newNR,I) == N) Neib(newNR,I) = E
			end do

			!------------------------ Deleting N -----------------------
			!print*,'Element: ',N,' Eliminated!',' : ',E

			do I=1,4
				Corn(N,I) = Corn(N,1)
			end do

			do I=1,4
				if(Corn(E,I) == C1) then
					Corn(E,I) = C2
					exit
				endif
			end do

			C1 = C2

		end do
		
	endif

end do
!Part 4:
if(done) then

	do I=1,4

		Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E,I),E)

	end do

endif

!===========================================================================================
End Subroutine NodeElimination
!*********************************************************************************************
