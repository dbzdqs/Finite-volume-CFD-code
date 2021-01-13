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
Subroutine ThreeEdgedNodesDividedByFourEdgedNode(Dim,Corn,Neib,X,Y,NBE,BFP,E1,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,E1
Intent(InOut)::Corn,Neib,X,Y,done

Integer::Dim,NBE,PC,E1,E2,E3,E4,EL,ER,TEC,QEC,SEC,Elm,A,Ai1,B,Bi1,Bi4,C,Ci2,D,Di2,Dir,E,Ei1,Ei2,F,Fi1,Fil,G,Gi,H,Hil,Hi4,I,Ii3,Ii4,J,Ji3,Jir,K,Ki,Na,Nb,Nc,Nd,Ne,Nf,Ng,Nh,Ni,Nj,Nk,P1,Pi1,P2,Pi2,Np1,Np2
Integer::L,M,N,getNeibour,Nde,Nef,Nfg,Ngh,Nhi,Nij,Njk,Nkd
Integer,Dimension(1:1000)::TElms,QElmsB,QElms,Points,SElms 
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent,SegmentJointIsPossible,done,IsInTheList,IsOnExteriorBoundary,A_Found,MoveCorners
Real(8)::x_value,y_value,xf,yf,xh,yh,xe,ye,xi,yi
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"
!--->>> Notice: Implementation of a series of Algorithms to remove THREE-EDGED NODES DIVIDED BY ONE FOUR-EDGED NODE. <<<---

do L=1,4 !------------------------------ Loop to find 'B' corner ---------------------------

	PC = 0
	A_Found = .False.
	done = .False.
    SEC = 0

	Bi1 = L !----------------------- Index of 'B' corner in E1 -----------------------------
	B = Corn(E1,L)

    if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,B)) then !--- Checking B on exterior boundary ---
!Part 1:    
	    Call GetSurroundingElements(Dim,Corn,Neib,B,E1,TElms,QElmsB,TEC,Nb)

	    if(Nb == 4) then !---------------- Condition of finding 'B' corner ---------------------

		    if(Bi1 == 1) then !-------- Finding Previous and Next corner to 'B' in E1 ----------

			    Pi1 = 2           !----------- Index of Next corner after B in E1 --------------
			    Pi2 = 4 		  !--------- Index of Previous corner after B in E1 ------------

			    P1 = Corn(E1,2)   !----------------- Next corner after B in E1 -----------------
			    P2 = Corn(E1,4)   !------------- Previous corner before B in E1 ----------------

		    elseif(Bi1 == 4) then

			    Pi1 = 1           !----------- Index of Next corner after B in E1 --------------
			    Pi2 = 3 		  !--------- Index of Previous corner after B in E1 ------------

			    P1 = Corn(E1,1)   !---------------- Next corner after B in E1 ------------------
			    P2 = Corn(E1,3)   !------------ Previous corner before B in E1 -----------------

		    else

			    Pi1 = Bi1+1       !----------- Index of Next corner after B in E1 --------------
			    Pi2 = Bi1-1 	  !--------- Index of Previous corner after B in E1 ------------

			    P1 = Corn(E1,Pi1) !---------------- Next corner after B in E1 ------------------
			    P2 = Corn(E1,Pi2) !------------ Previous corner before B in E1 -----------------

		    endif !----------- End of Finding Previous and Next corner to 'B' in E1 ------------

	       !------------------------- Finding Valances of P1 and P2 ----------------------------

		    Call GetSurroundingElements(Dim,Corn,Neib,P1,E1,TElms,QElms,TEC,Np1)
		    Call GetSurroundingElements(Dim,Corn,Neib,P2,E1,TElms,QElms,TEC,Np2)

	       if(Np1 == 3 .And. Np2 > 3) then !-------- Finding 'A' and 'E' corners in E1 ---------

			    Ai1 = Pi1 !------------------ Index of 'A' corner in E1  -----------------------
			    A = Corn(E1,Ai1)
			    A_Found = .True.

			    Ei1 = Pi2        !---------------- Index of 'E' corner in E1 -------------------
			    E = Corn(E1,Ei1) !---------- 'E' corner (opposite corner to 'A' in E1) ---------
			    Ne = Np2

			    if(.Not. IsInTheList(Points,PC,E)) then

				    PC = PC + 1
				    Points(PC) = E

			    endif

	       elseif(Np2 == 3 .And. Np1 > 3) then

			    Ai1 = Pi2 !------------------ Index of 'A' corner in E1 ------------------------
			    A = Corn(E1,Ai1)
			    A_Found = .True.

			    Ei1 = Pi1        !---------------- Index of 'E' corner in E1 -------------------
			    E = Corn(E1,Ei1) !---------- 'E' corner (opposite corner to 'A' in E1) ---------
			    Ne = Np1

			    if(.Not. IsInTheList(Points,PC,E)) then

				    PC = PC + 1
				    Points(PC) = E

			    endif

	       endif !----------------- End of Finding 'A' and 'E' corners in E1 -------------------

	       if(A_Found) then !------------------------- A is accepted ---------------------------

		       !------------------- Finding 'F' corner in E1 and its Valance -----------------------

		       do M=1,4
		   
				    if(Corn(E1,M)/=A .And. Corn(E1,M)/=B .And. Corn(E1,M)/=E) then

					    Fi1 = M !-------------------- Index of F in E1 -----------------------------
					    F = Corn(E1,M)

					    if(.Not. IsInTheList(Points,PC,F)) then

						    PC = PC + 1
						    Points(PC) = F

					    endif

					    exit

				    endif

		       end do

		       Call GetSurroundingElements(Dim,Corn,Neib,F,E1,TElms,QElms,TEC,Nf)

		       if(Nf > 3) then !---------------- Necessary condition on F valance ------------------

				    if(Ei1 > Bi1) then !------------------ Finding E2 Element ----------------------

					    if(Ei1 - Bi1 == 1) then

						    E2 = Neib(E1,Bi1)

					    else

						    E2 = Neib(E1,Ei1)

					    endif

				    else

					    if(Bi1 - Ei1 == 1) then

						    E2 = Neib(E1,Ei1)

					    else

						    E2 = Neib(E1,Bi1)

					    endif

				    endif !---------------------- End of Finding E2 Element ------------------------

				    !--------------------- Finding and Testing 'C' valance -------------------------

				    do M=1,4

					    if(areAdjacent(Dim,Corn,B,Corn(E2,M),E2) .And. Corn(E2,M) /= E) then

						    Ci2 = M !------------------ Index of 'C' in E2 -------------------------
						    C = Corn(E2,Ci2)
						    exit

					    endif

				    end do

				    Call GetSurroundingElements(Dim,Corn,Neib,C,E2,TElms,QElms,TEC,Nc)

				    if(Nc == 3) then !------------- Necessary condition on C valance ---------------

					    !------------------- Finding and Testing 'D' valance -----------------------

					    do M=1,4

						    if(areAdjacent(Dim,Corn,C,Corn(E2,M),E2) .And. Corn(E2,M) /= B) then

							    Di2 = M !---------------- Index of 'D' in E2 -----------------------
							    D = Corn(E2,Di2)
							    exit

						    endif

					    end do
					
					    Call GetSurroundingElements(Dim,Corn,Neib,D,E2,TElms,QElms,TEC,Nd)
					
					    if(Nd > 3) then !------------ Necessary condition on D valance -------------

						    !------------ Finding EL (Left Most Element in the pattern) ------------
						
						    if(Fi1 > Ai1) then

							    if(Fi1 - Ai1 == 1) then

								    EL = Neib(E1,Ai1)

							    else

								    EL = Neib(E1,Fi1)

							    endif
						
						    else

							    if(Ai1 - Fi1 == 1) then

								    EL = Neib(E1,Fi1)

							    else

								    EL = Neib(E1,Ai1)

							    endif
						
						    endif
						
						    !----------------- Finding 'G' and Testing its Valance -----------------

						    do M=1,4

							    if(areAdjacent(Dim,Corn,F,Corn(EL,M),EL) .And. Corn(EL,M) /= A) then

								    Gi = M !----------- Index of 'G' in EL ---------------
								    G = Corn(EL,Gi)

								    if(.Not. IsInTheList(Points,PC,G)) then

									    PC = PC + 1
									    Points(PC) = G

								    endif

								    exit

							    endif

						    end do
						
						    Call GetSurroundingElements(Dim,Corn,Neib,G,EL,TElms,QElms,TEC,Ng)
						
						    if(Ng > 3) then !------------ Necessary condition on G valance ---------

							    !--------------- Finding 'H' and Testing its Valance ---------------

							    do M=1,4

								    if(areAdjacent(Dim,Corn,G,Corn(EL,M),EL) .And. Corn(EL,M) /= F) then

									    Hil = M !------------- Index of 'H' in EL ------------------
									    H = Corn(EL,Hil)

									    if(.Not. IsInTheList(Points,PC,H)) then

										    PC = PC + 1
										    Points(PC) = H

									    endif

									    exit

								    endif

							    end do
							
							    Call GetSurroundingElements(Dim,Corn,Neib,H,EL,TElms,QElms,TEC,Nh)

							    if(Nh > 3) then !-------- Necessary condition on H valance ---------

								    !---------- Finding E4 (Adjacent Element to E1 and EL) ---------
						
								    if(Bi1 > Ai1) then

									    if(Bi1 - Ai1 == 1) then

										    E4 = Neib(E1,Ai1)

									    else

										    E4 = Neib(E1,Bi1)

									    endif
								
								    else

									    if(Ai1 - Bi1 == 1) then

										    E4 = Neib(E1,Bi1)

									    else

										    E4 = Neib(E1,Ai1)

									    endif
								
								    endif

								    !------------- Finding 'I' and Testing its Valance -------------

								    do M=1,4

									    if(areAdjacent(Dim,Corn,B,Corn(E4,M),E4) .And. Corn(E4,M) /= A) then

										    Ii4 = M !----------- Index of 'I' in E4 ----------------
										    I = Corn(E4,Ii4)

										    if(.Not. IsInTheList(Points,PC,I)) then

											    PC = PC + 1
											    Points(PC) = I

										    endif

										    exit

									    endif

								    end do
								
								    Call GetSurroundingElements(Dim,Corn,Neib,I,E4,TElms,QElms,TEC,Ni)

								    if(Ni > 3) then !------ Necessary condition on I valance -------

									    !-------- Finding E3 (Adjacent Element to E4 and E2) -------

									    do M=1,4

										    if(Corn(E4,M) == B) then

											    Bi4 = M
											    exit

										    endif

									    end do

									    if(Bi4 > Ii4) then

										    if(Bi4 - Ii4 == 1) then

											    E3 = Neib(E4,Ii4)

										    else

											    E3 = Neib(E4,Bi4)

										    endif

									    else

										    if(Ii4 - Bi4 == 1) then

											    E3 = Neib(E4,Bi4)

										    else

											    E3 = Neib(E4,Ii4)

										    endif

									    endif

									    !----------- Finding 'J' and Testing its Valance -----------

									    do M=1,4

										    if(areAdjacent(Dim,Corn,C,Corn(E3,M),E3) .And. Corn(E3,M) /= B) then

											    Ji3 = M !--------- Index of 'J' in E3 --------------
											    J = Corn(E3,Ji3)

											    if(.Not. IsInTheList(Points,PC,J)) then

												    PC = PC + 1
												    Points(PC) = J

											    endif

											    exit

										    endif

									    end do
									
									    Call GetSurroundingElements(Dim,Corn,Neib,J,E3,TElms,QElms,TEC,Nj)

									    if(Nj > 3) then !----- Necessary condition on J valance ----

										    !--- Finding ER (Right Most Element in the pattern) ----

										    if(Di2 > Ci2) then

											    if(Di2 - Ci2 == 1) then

												    ER = Neib(E2,Ci2)

											    else

												    ER = Neib(E2,Di2)

											    endif

										    else

											    if(Ci2 - Di2 == 1) then

												    ER = Neib(E2,Di2)

											    else

												    ER = Neib(E2,Ci2)

											    endif

										    endif
										
										    !----------- Finding 'K' and Testing its Valance -----------

										    do M=1,4

											    if(areAdjacent(Dim,Corn,D,Corn(ER,M),ER) .And. Corn(ER,M) /= C) then

												    Ki = M !---------- Index of 'K' in ER --------------
												    K = Corn(ER,Ki)

												    if(.Not. IsInTheList(Points,PC,K)) then

													    PC = PC + 1
													    Points(PC) = K

												    endif

												    exit

											    endif

										    end do
										
										    Call GetSurroundingElements(Dim,Corn,Neib,K,ER,TElms,QElms,TEC,Nk)
										
										    if(Nk > 3) then !------- Necessary condition on K valance ----------

											    !-------- Adding Last corner (D) of the RING of Elements -------

											    if(.Not. IsInTheList(Points,PC,D)) then

												    PC = PC + 1
												    Points(PC) = D

											    endif
!Part 2:
											    !======================== Specifying Neibour Elements =====================

											    Nde = getNeibour(Dim,Corn,Neib,D,E,E2)
                                                Call addToList(SElms,SEC,Nde)
											    Nef = getNeibour(Dim,Corn,Neib,E,F,E1)
                                                Call addToList(SElms,SEC,Nef)
											    Nfg = getNeibour(Dim,Corn,Neib,F,G,EL)
                                                Call addToList(SElms,SEC,Nfg)
											    Ngh = getNeibour(Dim,Corn,Neib,G,H,EL)
                                                Call addToList(SElms,SEC,Ngh)
											    Nhi = getNeibour(Dim,Corn,Neib,H,I,E4)
                                                Call addToList(SElms,SEC,Nhi)
											    Nij = getNeibour(Dim,Corn,Neib,I,J,E3)
                                                Call addToList(SElms,SEC,Nij)
											    Njk = getNeibour(Dim,Corn,Neib,J,K,ER)
                                                Call addToList(SElms,SEC,Njk)
											    Nkd = getNeibour(Dim,Corn,Neib,K,D,ER)
                                                Call addToList(SElms,SEC,Nkd)

											    !************>>>>>>>>>>>>>>> Implementation of 4 Different Cases <<<<<<<<<<<<<<<<<<************  

											    if((Nf >= 5 .And. Nj >= 5) .Or. (Nd >= 5 .And. Nh >= 5)) then !------>>> State1: Situation explained in Fig.9 in reference <<<------

												    if(Nf >= 5 .And. Nj >= 5) then
!Part 3:
													    if(Nd == 4 .And. Ne == 4 .And. Ng == 4 .And. Nh == 4 .And. Ni == 4 .And. Nk == 4) then

														    if(SegmentJointIsPossible(Dim,X,Y,E,H,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,D,I,Points,PC)) then

															    done = .True.

															    !------------------------------- Modifying EL --------------------------------

															    do M=1,4

																    if(Corn(EL,M) == A) then

																	    Corn(EL,M) = E !---------------- Mapping A to E ----------------------
																	    exit

																    endif

															    end do

															    !------------------------------- Modifying E2 --------------------------------

															    do M=1,4

																    if(Corn(E2,M) == B) then

																	    Corn(E2,M) = H !---------------- Mapping B to H ----------------------

																    endif

																    if(Corn(E2,M) == C) then

																	    Corn(E2,M) = I !---------------- Mapping C to I ----------------------

																    endif

															    end do

															    !------------------------------- Modifying ER --------------------------------

															    do M=1,4

																    if(Corn(ER,M) == C) then

																	    Corn(ER,M) = I !---------------- Mapping C to I ----------------------
																	    exit

																    endif

															    end do

															    !--------------------------- Modifying EL Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,EL,E,F,Nef)
															    Call setNeibour(Dim,Corn,Neib,EL,H,E,E2)

															    !--------------------------- Modifying E2 Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,E2,E,H,EL)
															    Call setNeibour(Dim,Corn,Neib,E2,H,I,Nhi)
															    Call setNeibour(Dim,Corn,Neib,E2,I,D,ER)

															    !--------------------------- Modifying ER Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,ER,D,I,E2)
															    Call setNeibour(Dim,Corn,Neib,ER,I,J,Nij)

															    !========================= Modifying Neibours of Neibours ====================

															    Call setNeibour(Dim,Corn,Neib,Nef,E,F,EL)
															    Call setNeibour(Dim,Corn,Neib,Nhi,H,I,E2)
															    Call setNeibour(Dim,Corn,Neib,Nij,I,J,ER)

															    !====================== Deleting E1,E3,E4 Elements ===========================

															    do M=2,4

																    Corn(E1,M) = Corn(E1,1)
																    Corn(E3,M) = Corn(E3,1)
																    Corn(E4,M) = Corn(E4,1)

															    end do

															    !====================== Smoothing EL,E2,ER Elements ==========================

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(EL,M),EL)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E2,M),E2)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ER,M),ER)

															    end do

															    exit

														    endif !--------------------- End of Segment Joint Checkings ----------------------

													    endif

												    else
!Part 4:
													    if(Ne == 4 .And. Nf == 4 .And. Ng == 4 .And. Ni == 4 .And. Nj == 4 .And. Nk == 4) then

														    if(SegmentJointIsPossible(Dim,X,Y,F,I,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,E,J,Points,PC)) then

															    done = .True.

															    !------------------------------- Modifying EL --------------------------------

															    do M=1,4

																    if(Corn(EL,M) == A) then

																	    Corn(EL,M) = I !---------------- Mapping A to I ----------------------
																	    exit

																    endif

															    end do

															    !------------------------------- Modifying E1 --------------------------------

															    do M=1,4

																    if(Corn(E1,M) == A) then

																	    Corn(E1,M) = I !---------------- Mapping A to I ----------------------

																    endif

																    if(Corn(E1,M) == B) then

																	    Corn(E1,M) = J !---------------- Mapping B to J ----------------------

																    endif

															    end do

															    !------------------------------- Modifying ER --------------------------------

															    do M=1,4

																    if(Corn(ER,M) == C) then

																	    Corn(ER,M) = E !---------------- Mapping C to E ----------------------
																	    exit

																    endif

															    end do

															    !--------------------------- Modifying EL Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,EL,H,I,Nhi)
															    Call setNeibour(Dim,Corn,Neib,EL,I,F,E1)

															    !--------------------------- Modifying E1 Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,E1,F,I,EL)
															    Call setNeibour(Dim,Corn,Neib,E1,I,J,Nij)
															    Call setNeibour(Dim,Corn,Neib,E1,J,E,ER)

															    !--------------------------- Modifying ER Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,ER,D,E,Nde)
															    Call setNeibour(Dim,Corn,Neib,ER,E,J,E1)

															    !========================= Modifying Neibours of Neibours ====================

															    Call setNeibour(Dim,Corn,Neib,Nde,D,E,ER)
															    Call setNeibour(Dim,Corn,Neib,Nhi,H,I,EL)
															    Call setNeibour(Dim,Corn,Neib,Nij,I,J,E1)

															    !====================== Deleting E2,E3,E4 Elements ===========================

															    do M=2,4

																    Corn(E2,M) = Corn(E2,1)
																    Corn(E3,M) = Corn(E3,1)
																    Corn(E4,M) = Corn(E4,1)

															    end do

															    !====================== Smoothing EL,E1,ER Elements ==========================

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(EL,M),EL)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E1,M),E1)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ER,M),ER)

															    end do

															    exit


														    endif !--------------------- End of Segment Joint Checkings ----------------------

													    endif

												    endif

											    elseif((Nd >= 5 .And. Nf >= 5 .And. Ni >= 5) .Or. (Ne >= 5 .And. Nh >= 5 .And. Nj >= 5)) then !---->>> State2: Situation explained in Fig.26 in reference <<<----

												    if(Nd >= 5 .And. Nf >= 5 .And. Ni >= 5) then

													    if(Ne == 4 .And. Ng == 4 .And. Nh == 4 .And. Nj == 4 .And. Nk == 4) then
!Part 5:
														    if(SegmentJointIsPossible(Dim,X,Y,E,H,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,E,J,Points,PC)) then

															    done = .True.

															    !------------------------------- Modifying EL --------------------------------

															    do M=1,4

																    if(Corn(EL,M) == A) then

																	    Corn(EL,M) = E !---------------- Mapping A to E ----------------------
																	    exit

																    endif

															    end do
															
															    !------------------------------- Modifying E3 --------------------------------

															    do M=1,4

																    if(Corn(E3,M) == B) then

																	    Corn(E3,M) = H !---------------- Mapping B to H ----------------------

																    endif

																    if(Corn(E3,M) == C) then

																	    Corn(E3,M) = E !---------------- Mapping C to E ----------------------

																    endif

															    end do
															
															    !------------------------------- Modifying ER --------------------------------

															    do M=1,4

																    if(Corn(ER,M) == C) then

																	    Corn(ER,M) = E !---------------- Mapping C to E ----------------------
																	    exit

																    endif

															    end do
															
															    !--------------------------- Modifying EL Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,EL,E,F,Nef)
															    Call setNeibour(Dim,Corn,Neib,EL,H,E,E3)
															
															    !--------------------------- Modifying E3 Neibours ---------------------------
															
															    Call setNeibour(Dim,Corn,Neib,E3,J,E,ER)
															    Call setNeibour(Dim,Corn,Neib,E3,E,H,EL)
															    Call setNeibour(Dim,Corn,Neib,E3,H,I,Nhi)
															
															    !--------------------------- Modifying ER Neibours ---------------------------
															
															    Call setNeibour(Dim,Corn,Neib,ER,D,E,Nde)
															    Call setNeibour(Dim,Corn,Neib,ER,E,J,E3)
															
															    !========================= Modifying Neibours of Neibours ====================

															    Call setNeibour(Dim,Corn,Neib,Nde,D,E,ER)
															    Call setNeibour(Dim,Corn,Neib,Nhi,H,I,E3)
															    Call setNeibour(Dim,Corn,Neib,Nef,E,F,EL)
															
															    !====================== Deleting E1,E2,E4 Elements ===========================

															    do M=2,4

																    Corn(E1,M) = Corn(E1,1)
																    Corn(E2,M) = Corn(E2,1)
																    Corn(E4,M) = Corn(E4,1)

															    end do
															
															    !====================== Smoothing EL,E3,ER Elements ==========================

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(EL,M),EL)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E3,M),E3)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ER,M),ER)

															    end do

															    exit	

														    endif !--------------------- End of Segment Joint Checkings ----------------------

													    endif

												    else 

													    if(Nd == 4 .And. Nf == 4 .And. Ng == 4 .And. Ni == 4 .And. Nk == 4) then
!Part 6:
														    if(SegmentJointIsPossible(Dim,X,Y,D,I,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,F,I,Points,PC)) then

															    done = .True.

															    !--------------------- Modifying EL (AFGH to IFGH) ---------------------------

															    do M=1,4

																    if(Corn(EL,M) == A) then

																	    Corn(EL,M) = I !---------------- Mapping A to I ----------------------
																	    exit

																    endif

															    end do

															    !----------------------- Modifying E1 (ABEF to IDEF) -------------------------

															    do M=1,4

																    if(Corn(E1,M) == B) then

																	    Corn(E1,M) = D !---------------- Mapping B to D ----------------------

																    endif

																    if(Corn(E1,M) == A) then

																	    Corn(E1,M) = I !---------------- Mapping A to I ----------------------

																    endif

															    end do

															    !--------------------- Modifying ER (DCJK to DIJK) ---------------------------

															    do M=1,4

																    if(Corn(ER,M) == C) then

																	    Corn(ER,M) = I !---------------- Mapping C to I ----------------------
																	    exit

																    endif

															    end do

															    !--------------------------- Modifying EL Neibours ---------------------------

															    Call setNeibour(Dim,Corn,Neib,EL,I,F,E1)
															    Call setNeibour(Dim,Corn,Neib,EL,H,I,Nhi)

															    !--------------------------- Modifying E1 Neibours ---------------------------
															
															    Call setNeibour(Dim,Corn,Neib,E1,D,E,Nde)
															    Call setNeibour(Dim,Corn,Neib,E1,F,I,EL)
															    Call setNeibour(Dim,Corn,Neib,E1,I,D,ER)
															
															    !--------------------------- Modifying ER Neibours --------------------------- 
															
															    Call setNeibour(Dim,Corn,Neib,ER,D,I,E1)
															    Call setNeibour(Dim,Corn,Neib,ER,I,J,Nij)

															    !========================= Modifying Neibours of Neibours ====================

															    Call setNeibour(Dim,Corn,Neib,Nde,D,E,E1)
															    Call setNeibour(Dim,Corn,Neib,Nhi,H,I,EL)
															    Call setNeibour(Dim,Corn,Neib,Nij,I,J,ER)
															
															    !====================== Deleting E2,E3,E4 Elements ===========================

															    do M=2,4

																    Corn(E2,M) = Corn(E2,1)
																    Corn(E3,M) = Corn(E3,1)
																    Corn(E4,M) = Corn(E4,1)

															    end do

															    !====================== Smoothing EL,E1,ER Elements ==========================

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(EL,M),EL)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E1,M),E1)

															    end do

															    do M=1,4

																    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ER,M),ER)

															    end do

															    exit

														    endif !--------------------- End of Segment Joint Checkings ----------------------

													    endif 

												    endif

											    elseif(Ne >= 5 .And. Ni >= 5) then !---->>> State3: Situation explained in Fig.27 in reference <<<----

												    if(Nd == 4 .And. Nf == 4 .And. Ng == 4 .And. Nh == 4 .And. Nj == 4 .And. Nk == 4) then
!Part 7:
													    if(SegmentJointIsPossible(Dim,X,Y,H,D,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,F,J,Points,PC)) then

														    done = .True.

														    !--------------------- Modifying EL (AFGH to BFGH) ---------------------------

														    do M=1,4

															    if(Corn(EL,M) == A) then

																    Corn(EL,M) = B !---------------- Mapping A to B ----------------------
																    exit

															    endif

														    end do

														    !--------------------- Modifying E1 (ABEF to BDEF) ---------------------------

														    do M=1,4

															    if(Corn(E1,M) == B) then

																    Corn(E1,M) = D !---------------- Mapping B to D ----------------------

															    endif

															    if(Corn(E1,M) == A) then

																    Corn(E1,M) = B !---------------- Mapping A to B ----------------------

															    endif

														    end do

														    !--------------------- Modifying E4 (AHIB to AHIJ) ---------------------------

														    do M=1,4

															    if(Corn(E4,M) == B) then

																    Corn(E4,M) = J !---------------- Mapping B to J ----------------------

															    endif

															    if(Corn(E4,M) == A) then

																    Corn(E4,M) = B !---------------- Mapping A to B ----------------------

															    endif

														    end do

														    !--------------------- Modifying ER (DCJK to DAJK) ---------------------------

														    do M=1,4

															    if(Corn(ER,M) == C) then

																    Corn(ER,M) = B !---------------- Mapping C to B ----------------------
																    exit

															    endif

														    end do

														    !--------------------------- Modifying E1 Neibours ---------------------------

														    Call setNeibour(Dim,Corn,Neib,E1,B,D,ER)
														    Call setNeibour(Dim,Corn,Neib,E1,D,E,Nde)

														    !--------------------------- Modifying E4 Neibours ---------------------------

														    Call setNeibour(Dim,Corn,Neib,E4,J,B,ER)
														    Call setNeibour(Dim,Corn,Neib,E4,I,J,Nij)

														    !--------------------------- Modifying ER Neibours ---------------------------

														    Call setNeibour(Dim,Corn,Neib,ER,D,B,E1)
														    Call setNeibour(Dim,Corn,Neib,ER,B,J,E4)

														    !========================= Modifying Neibours of Neibours ====================

														    Call setNeibour(Dim,Corn,Neib,Nde,D,E,E1)
														    Call setNeibour(Dim,Corn,Neib,Nij,I,J,E4)

														    !========================= Deleting E2,E3 Elements ===========================

														    do M=2,4

															    Corn(E2,M) = Corn(E2,1)
															    Corn(E3,M) = Corn(E3,1)

														    end do

														    !===================== Smoothing EL,E1,E4,ER Elements ========================

														    do M=1,4

															    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(EL,M),EL)

														    end do

														    do M=1,4

															    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E1,M),E1)

														    end do

														    do M=1,4

															    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(E4,M),E4)

														    end do

														    do M=1,4

															    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ER,M),ER)

														    end do

														    exit

													    endif !--------------------- End of Segment Joint Checkings ----------------------

												    endif

											    elseif(Ng >= 5 .Or. Nk >= 5) then !---->>> State4: Situation explained in Fig.28 in reference <<<----

												    if(Nd == 4 .And. Ne == 4 .And. Nf == 4 .And. Nh == 4 .And. Ni == 4 .And. Nj == 4) then
!Part 8:
													    if(SegmentJointIsPossible(Dim,X,Y,F,H,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,E,I,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,D,J,Points,PC)) then

														    done = .True.
                                                            
                                                            xf = X(F)
                                                            yf = Y(F)
                                                            
                                                            xh = X(H)
                                                            yh = Y(H)

														    if(MoveCorners(Dim,NBE,BFP,Corn,X,Y,F,H,SElms,SEC)) then !------------- Moving F to midpoint of FH -------------
                                                                
                                                                xe = X(E)
                                                                ye = Y(E)
                                                                
                                                                xi = X(I)
                                                                yi = Y(I)
                                                                
                                                                if(MoveCorners(Dim,NBE,BFP,Corn,X,Y,E,I,SElms,SEC)) then !------ Moving E to midpoint of EI ---------
                                                                
                                                                    if(MoveCorners(Dim,NBE,BFP,Corn,X,Y,D,J,SElms,SEC)) then !-------- Moving D to midpoint of DJ ---------
                                                                   
                                                                        !------------------- Mapping H to F -------------------
														
														                Call GetSurroundingElements(Dim,Corn,Neib,H,EL,TElms,QElms,TEC,QEC)
														
														                do M=1,QEC
														
															                Elm = QElms(M)

															                do N=1,4

																                if(Corn(Elm,N) == H) then

																	                Corn(Elm,N) = F
																	                exit

																                endif

															                end do

														                end do
														
														                !------------------- Mapping I to E -------------------
														
														                Call GetSurroundingElements(Dim,Corn,Neib,I,E4,TElms,QElms,TEC,QEC)
														
														                do M=1,QEC
														
															                Elm = QElms(M)

															                do N=1,4

																                if(Corn(Elm,N) == I) then

																	                Corn(Elm,N) = E
																	                exit

																                endif

															                end do

														                end do 

														                !------------------- Mapping J to D -------------------
														
														                Call GetSurroundingElements(Dim,Corn,Neib,J,ER,TElms,QElms,TEC,QEC)
														
														                do M=1,QEC
														
															                Elm = QElms(M)

															                do N=1,4

																                if(Corn(Elm,N) == J) then

																	                Corn(Elm,N) = D
																	                exit

																                endif

															                end do

														                end do

														                !------------- Modifying Neibours of Neibours ---------

														                Call setNeibour(Dim,Corn,Neib,Nkd,D,K,Njk)
														                Call setNeibour(Dim,Corn,Neib,Njk,D,K,Nkd)

														                Call setNeibour(Dim,Corn,Neib,Nde,D,E,Nij)
														                Call setNeibour(Dim,Corn,Neib,Nij,D,E,Nde)

														                Call setNeibour(Dim,Corn,Neib,Nef,E,F,Nhi)
														                Call setNeibour(Dim,Corn,Neib,Nhi,E,F,Nef)

														                Call setNeibour(Dim,Corn,Neib,Nfg,F,G,Ngh)
														                Call setNeibour(Dim,Corn,Neib,Ngh,F,G,Nfg)

														                !------------------ Deleting Elements -----------------

														                do M=2,4
															
															                Corn(EL,M) = Corn(EL,1)
															                Corn(E1,M) = Corn(E1,1)
															                Corn(E2,M) = Corn(E2,1)
															                Corn(E3,M) = Corn(E3,1)
															                Corn(E4,M) = Corn(E4,1)
															                Corn(ER,M) = Corn(ER,1)

														                end do

														                !===================== Smoothing Neibour Elements =======================

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nde,M),Nde)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nef,M),Nef)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nfg,M),Nfg)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Ngh,M),Ngh)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nhi,M),Nhi)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nij,M),Nij)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Njk,M),Njk)

														                end do

														                do M=1,4

															                Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nkd,M),Nkd)

														                end do

														                exit
                                                                        
                                                                    else
                                                                        
                                                                        X(F) = xf
                                                                        Y(F) = yf
                                                            
                                                                        X(H) = xh
                                                                        Y(H) = yh
                                                                        
                                                                        X(E) = xe
                                                                        Y(E) = ye
                                                                
                                                                        X(I) = xi
                                                                        Y(I) = yi
                                                                        
                                                                        done = .False.
                                                                        
                                                                    endif
                                                                    
                                                                else
                                                                
                                                                    X(F) = xf
                                                                    Y(F) = yf
                                                            
                                                                    X(H) = xh
                                                                    Y(H) = yh
                                                                    
                                                                    done = .False.
                                                                    
                                                                endif
                                                                
                                                            else
                                                                
                                                                done = .False.

                                                            endif

													    endif !--------------------- End of Segment Joint Checkings ----------------------

												    endif

											    endif !*************>>>>>>>>>>>>> End of Implementation of Cases <<<<<<<<<<<<<<<<<************  
										
										    endif !-- End of Necessary condition on K valance ------	

									    endif !------ End of Necessary condition on J valance ------

								    endif !---------- End of Necessary condition on I valance ------

							    endif !-------------- End of Necessary condition on H valance ------
						
						    endif !------------------ End of Necessary condition on G valance ------		
					
					    endif !---------------------- End of Necessary condition on D valance ------	

				    endif !-------------------------- End of Necessary condition on C valance ------

		       endif !------------------------------- End of Necessary condition on F valance ------

	       endif !----------------------------------- End of A accept ------------------------------

        endif !-------------------------------------- End of Condition of finding 'B' corner -------
    
    endif !-------------------------------------- End of Checking B on the exterior boundary ----

end do !----------------------------------------- End of Main Loop -----------------------------

!===========================================================================================
End Subroutine ThreeEdgedNodesDividedByFourEdgedNode
!*********************************************************************************************
