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
!// Date: Nov., 15, 2014                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine ThreeFiveOppoThree(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,TEC,ME,I,J,K,A,Ai,B,Bi,C,Ci,D,Di,Na,Nb,Nc1,Nd,P1,Np1,index,Elm,newNode,newElement,getNeibour,Nab,Nbc,Cnt
Integer,Dimension(1:1000)::QElmsA,TElms,QElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::isQuadNode,Pattern_Found,done,done1,IsDiagonalInsideQuad,areAdjacent
Real(8)::QuadQuality
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do I=1,4

	Pattern_Found = .False.
	done = .False.

	Ai = I
	A = Corn(ME,Ai)
!Part 1:
	Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElmsA,TEC,Na)

	if(Na == 5) then !------------------------ Condition on A Valance ----------------------

		Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

		Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc1)

		if(Nc1 == 3) then !--------------------- Condition on C Valance --------------------

			!---------------------- Finding Next and Previous Corner to 'A' ----------------

			Call getNextCorner(Dim,Corn,ME,A,B,Bi)
			Call getPrevCorner(Dim,Corn,ME,A,D,Di)

			Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElms,TEC,Nb)
			Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)

			if(Nb == 4 .And. Nd == 4) then !------- Checking 'B' and 'D' Valances ----------

				!-- Finding an Element among 'C' Elements having an adjacent corner to 'C' with LOW Valance --
				
				P1 = B
				Elm = ME
				Cnt = 0

				do

					do K=1,4

						if(areAdjacent(Dim,Corn,A,Corn(Elm,K),Elm) .And. Corn(Elm,K) /= P1) then

							Cnt = Cnt + 1
							 	
							P1 = Corn(Elm,K)

							if(Cnt /= 3) then

								Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)

							endif

							exit

						endif

					end do

					if(Cnt == 3) exit

				end do
				
				Call GetSurroundingElements(Dim,Corn,Neib,P1,Elm,TElms,QElms,TEC,Np1)
				
				if(Np1 == 3) then
				
					Pattern_Found = .True.	

				endif	

				if(Pattern_Found) then
!Part 2:
					!------- Adding new Node to Center of ME (AC Must be inside ME) --------

					if(IsDiagonalInsideQuad(Dim,Corn,X,Y,ME,A,C) .And. (QuadQuality(Dim,X,Y,A,B,C,D)/=0)) then

						NP = NP + 1
						newNode = NP

						X(newNode) = (X(A) + X(C))/2
						Y(newNode) = (Y(A) + Y(C))/2

						!--------------- Finding Nab and Nbc Elements (ME naibours) --------

						Nab = getNeibour(Dim,Corn,Neib,A,B,ME)
						Nbc = getNeibour(Dim,Corn,Neib,B,C,ME)

						!----------------------- Adding newElement -------------------------

						NC = NC + 1
						newElement = NC

						Corn(newElement,Ai) = A
						Corn(newElement,Bi) = B
						Corn(newElement,Ci) = C
						Corn(newElement,Di) = newNode

						!--------------------------- Modifying ME --------------------------

						Corn(ME,Ai) = A
						Corn(ME,Bi) = newNode
						Corn(ME,Ci) = C
						Corn(ME,Di) = D

						!------------------- Specifying newElement Neibours ----------------

						Call setNeibour(Dim,Corn,Neib,newElement,A,B,Nab)
						Call setNeibour(Dim,Corn,Neib,newElement,B,C,Nbc)
						Call setNeibour(Dim,Corn,Neib,newElement,newNode,C,ME)
						Call setNeibour(Dim,Corn,Neib,newElement,newNode,A,ME)

						!------------------ Modifying Neibours of newElement ---------------

						Call setNeibour(Dim,Corn,Neib,Nab,A,B,newElement)
						Call setNeibour(Dim,Corn,Neib,Nbc,B,C,newElement)
						Call setNeibour(Dim,Corn,Neib,ME,newNode,A,newElement)
						Call setNeibour(Dim,Corn,Neib,ME,newNode,C,newElement)
!Part 3:
						!------------------------ Resolving Pattern ------------------------

						Call ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,newNode,A,P1,ME,done)
						
						if(done) then
						
							exit

						else
!Part 4:						
							Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,ME,done1)
							
							NC = NC - 1

						endif 

					endif !-------------------- End of AC necessary conditions ------------- 

				endif

			endif !----------------- End of Checking 'B' and 'D' Valances ------------------

		endif !------------------------ End of Condition on C Valance ---------------------- 

	endif !------------------------- End of Condition on A Valance -------------------------

end do !-------------------------- End of Main Loop ('I' Loop) -----------------------------

!===========================================================================================
End Subroutine ThreeFiveOppoThree
!*********************************************************************************************
