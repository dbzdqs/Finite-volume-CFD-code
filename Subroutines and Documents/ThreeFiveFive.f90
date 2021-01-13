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
Subroutine ThreeFiveFive(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,ME,done) !--->>> Pattern in Fig 33 of reference
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,TEC,ME,NE1,NE2,N1,N2,index,newElement,newP,newE,I,J,A,Ai,B,Bi,C,Ci,D,Di,E,F,P1,P2,P3,P4,P5,Np3,Np4,Np5,Na,Nb,Nc1,Nd,Ne,Nf,Elm,getNeibour
Integer,Dimension(1:1000)::QElms,TElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsOnExteriorBoundary,done,done1,done2,areAdjacent,IsDiagonalInsideQuad,ElementFound
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do I=1,4

	done = .False.
	done1 = .False.
	done2 = .False.
	ElementFound = .False.
!Part 1:
	Ai = I
	A = Corn(ME,Ai)

	if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,A)) then !---- End of A Boundary Checking ----
	
		Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElms,TEC,Na)

		if(Na == 5) then !---------------------- Condition on A Valance --------------------

			Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

			if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,C)) then !- End of C Boundary Checking -

				Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc1)

				if(Nc1 == 4) then !----------------- Condition on C Valance ----------------

					Call getNextCorner(Dim,Corn,ME,A,B,Bi)
					Call getPrevCorner(Dim,Corn,ME,A,D,Di)

					if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,B) .And. .Not. IsOnExteriorBoundary(Dim,NBE,BFP,D)) then

						Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElms,TEC,Nb)
						Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)

						if(Nb == 3 .And. Nd == 4) then

							ElementFound = .True.

							P1 = B
							P2 = D

						elseif(Nb == 4 .And. Nd == 3) then

							ElementFound = .True.

							P1 = D
							P2 = B

						endif

					endif

					if(ElementFound) then !--------- Condition on ElementFound ------------- 

						NE1 = getNeibour(Dim,Corn,Neib,A,P2,ME)

						do J=1,4

							if(areAdjacent(Dim,Corn,A,Corn(NE1,J),NE1) .And. Corn(NE1,J) /= P2) then
							
								P3 = Corn(NE1,J)
								exit 

							endif	

						end do

						if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,P3)) then !-- P3 Boundary Check --

							Call GetSurroundingElements(Dim,Corn,Neib,P3,NE1,TElms,QElms,TEC,Np3)

							if(Np3 == 5) then !-------- Condition on P3 Valence ------------

								Call getOppoCorner(Dim,Corn,NE1,A,P4,index)

								if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,P4)) then !-- P4 Boundary Check --
								
									Call GetSurroundingElements(Dim,Corn,Neib,P4,NE1,TElms,QElms,TEC,Np4)
									
									if(Np4 == 4) then !----- Condition on P4 Valence -------
									
										NE2 = getNeibour(Dim,Corn,Neib,P3,P4,NE1)

										do J=1,4

											if(areAdjacent(Dim,Corn,P3,Corn(NE2,J),NE2) .And. Corn(NE2,J) /= P4) then
											
												P5 = Corn(NE2,J)
												exit 

											endif	

										end do

										if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,P4)) then !-- P4 Boundary Check --

											Call GetSurroundingElements(Dim,Corn,Neib,P5,NE2,TElms,QElms,TEC,Np5)

											if(Np5 <= 4) then !------- Pattern Found -------
!Part 2:
												if(IsDiagonalInsideQuad(Dim,Corn,X,Y,NE1,P2,P3)) then

													NP = NP + 1
													newP = NP

													X(newP) = (X(P2) + X(P3))/2
													Y(newP) = (Y(P2) + Y(P3))/2

													N1 = getNeibour(Dim,Corn,Neib,P2,P4,NE1)
													N2 = getNeibour(Dim,Corn,Neib,A,P3,NE1)

													!---------- Adding new Element ---------

													NC = NC + 1
													newE = NC

													do J=1,4

														if(Corn(NE1,J) == A) Corn(newE,J) = A
														if(Corn(NE1,J) == P2) Corn(newE,J) = P2
														if(Corn(NE1,J) == P3) Corn(newE,J) = P3
														if(Corn(NE1,J) == P4) Corn(newE,J) = newP

													end do

													!------------- Modifying NE1 -----------

													do J=1,4

														if(Corn(NE1,J) == A) then

															Corn(NE1,J) = newP
															exit

														endif

													end do

													!-------- Specifying newE Neibours -----

													Call setNeibour(Dim,Corn,Neib,newE,A,P2,ME)
													Call setNeibour(Dim,Corn,Neib,newE,A,P3,N2)
													Call setNeibour(Dim,Corn,Neib,newE,P2,newP,NE1)
													Call setNeibour(Dim,Corn,Neib,newE,P3,newP,NE1)

													!-------- Specifying NE1 Neibours ------

													Call setNeibour(Dim,Corn,Neib,NE1,P2,newP,newE)
													Call setNeibour(Dim,Corn,Neib,NE1,P3,newP,newE)

													!----- Specifying Neibours Neibours ----

													Call setNeibour(Dim,Corn,Neib,N2,A,P3,newE)
													Call setNeibour(Dim,Corn,Neib,ME,A,P2,newE)
!Part 3:
													Call CollapseOperation(Dim,Corn,Neib,X,Y,NBE,BFP,P1,A,P2,ME,done) !--- Notice: Node P1 is Removed

													if(done) then
!Part 4:
														Call ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,newP,P2,C,NE1,done1)
														Call ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,newP,P3,P5,NE1,done2)

														if(.Not. done1 .And. .Not. done2) then

															Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,NE1,done1)
												
															NC = NC - 1

														else

															done = .True.
															exit

														endif

													else

														Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,NE1,done1)
												
														NC = NC - 1

													endif

												endif

											endif !--------- End of Pattern Found ----------

										endif !-------- End of P5 Boundary Checking --------

									endif !-------- End of Condition on P4 Valence ---------
								
								endif !----------- End of P4 Boundary Checking -------------	

							endif !------------ End of Condition on P3 Valence -------------

						endif !-------------- End of P3 Boundary Checking ------------------ 

					endif !-------------- End of Condition on ElementFound -----------------

				endif !------------------- End of Condition on C Valance -------------------

			endif !----------------------- End of C Boundary Checking ----------------------

		endif !------------------------ End of Condition on A Valance ----------------------

	endif !--------------------------- End of A Boundary Checking --------------------------

end  do !---------------------------- End of Main Loop (I Loop) ----------------------------

!===========================================================================================
End Subroutine ThreeFiveFive
!*********************************************************************************************
