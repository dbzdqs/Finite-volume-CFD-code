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
Subroutine SegmentElimination(Dim,Corn,Neib,X,Y,NBE,BFP,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::Corn,Neib,X,Y,done

Integer::Dim,NBE,ME,NE,NL,NR,I,J,K,PC,Na,Nc,Nd,Ne1,Nf,Ng,Nh,P1,P2,A,B,C,D,E,F,G,H,TEC,QEC,Elm,Ai1,Ai2,Bi1,Bi2,Ci,Di,Ei,Fi,Gi,Hi,N1,N2,Pi1,Pi2,getNeibour
Integer::index1,index2,Ncd,Nde,Nef,Nfg,Ngh,Nhc,CalcNodeValance,val
Integer,Dimension(1:1000)::TElms,QElms,Points
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent,SegmentJointIsPossible,done,IsInTheList,IsOnExteriorBoundary,SegmentFound
Real(8)::Q1,Q2,QuadQuality,x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: 1. Lee, C.K.; Lo, S.H. (1994) A new scheme for the generation of a graded quadrilateral mesh.
!--->>> Reference: 2. Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

done = .False.

do I=1,4
!Part 1:	
	SegmentFound = .False.

	PC = 0

	Ai1 = I
	A = Corn(ME,Ai1)

	if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,A)) then !- Checking A not to be on the exterior boundary 

		Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElms,TEC,Na)

		if(Na == 3) then !---------------------- Initial Condition -----------------------------

			if(Ai1 == 1) then

				Pi1 = 2         !--------------- Index of Next corner after A ------------------
				Pi2 = 4 		!------------- Index of Previous corner after A ----------------

				P1 = Corn(ME,2) !-------------------- Next corner after A ----------------------
				P2 = Corn(ME,4) !---------------- Previous corner before A ---------------------

			elseif(Ai1 == 4) then

				Pi1 = 1         !--------------- Index of Next corner after A ------------------
				Pi2 = 3 		!------------- Index of Previous corner after A ----------------

				P1 = Corn(ME,1) !-------------------- Next corner after A ----------------------
				P2 = Corn(ME,3) !---------------- Previous corner before A ---------------------

			else

				Pi1 = Ai1+1     !--------------- Index of Next corner after A ------------------
				Pi2 = Ai1-1 	!------------- Index of Previous corner after A ----------------

				P1 = Corn(ME,Ai1+1) !------------------ Next corner after A --------------------
				P2 = Corn(ME,Ai1-1) !-------------- Previous corner before A -------------------

			endif

			Call GetSurroundingElements(Dim,Corn,Neib,P1,ME,TElms,QElms,TEC,N1)
			Call GetSurroundingElements(Dim,Corn,Neib,P2,ME,TElms,QElms,TEC,N2)

			if(N1 == 3 .And. N2 /= 3) then !--- Next node to A satisfies neccessary condition --
			
				if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,P1)) then

					SegmentFound = .True.
					Bi1 = Pi1
					B = P1

				endif	

			elseif(N2 == 3 .And. N1 /= 3) then !- Prev node to A satisfies neccessary condition -

				if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,P2)) then

					SegmentFound = .True.
					Bi1 = Pi2
					B = P2

				endif

			endif
 
			if(SegmentFound) then !------------------ Segment Found condition ----------------------
!Part 2:
				!---------------- Finding NE (Neibour of ME having AB segment in common) -----------
				
                NE = getNeibour(Dim,Corn,Neib,A,B,ME)

				!----------------------------- Finding A and B indices in NE -----------------------

				do J=1,4
					
					if(Corn(NE,J) == A) Ai2 = J
					if(Corn(NE,J) == B) Bi2 = J 

				end do

				!------------------------------ Finding other ME corners ---------------------------

				do J=1,4

					if(areAdjacent(Dim,Corn,B,Corn(ME,J),ME) .And. Corn(ME,J) /= A) then

						Ci = J
						C = Corn(ME,J)

						if(.Not. IsInTheList(Points,PC,C)) then

							PC = PC + 1
							Points(PC) = C

						endif

						exit

					endif

				end do

				do J=1,4

					if(areAdjacent(Dim,Corn,A,Corn(ME,J),ME) .And. Corn(ME,J) /= B) then

						Di = J
						D = Corn(ME,J)

						if(.Not. IsInTheList(Points,PC,D)) then

							PC = PC + 1
							Points(PC) = D

						endif

						exit

					endif

                end do

				!------------ Finding NL and NR (common Quad neibours of both ME and NE) -----------

                NL = getNeibour(Dim,Corn,Neib,A,D,ME)
                NR = getNeibour(Dim,Corn,Neib,B,C,ME)

				!------------------ Finding E in NL (opposite corner to A in NL) -------------------

				do J=1,4

					if(areAdjacent(Dim,Corn,D,Corn(NL,J),NL) .And. Corn(NL,J) /= A) then

						E = Corn(NL,J)

						if(.Not. IsInTheList(Points,PC,E)) then

							PC = PC + 1
							Points(PC) = E

						endif

						exit

					endif

				end do

				!------------------------------ Finding other NE corners ---------------------------

				do J=1,4

					if(areAdjacent(Dim,Corn,A,Corn(NE,J),NE) .And. Corn(NE,J) /= B) then

						Fi = J
						F = Corn(NE,J)

						if(.Not. IsInTheList(Points,PC,F)) then

							PC = PC + 1
							Points(PC) = F

						endif

						exit

					endif

				end do

				do J=1,4

					if(areAdjacent(Dim,Corn,B,Corn(NE,J),NE) .And. Corn(NE,J) /= A) then

						Gi = J
						G = Corn(NE,J)

						if(.Not. IsInTheList(Points,PC,G)) then

							PC = PC + 1
							Points(PC) = G

						endif

						exit

					endif

				end do

				!------------------ Finding H in NR (opposite corner to B in NR) -------------------

				do J=1,4

					if(areAdjacent(Dim,Corn,C,Corn(NR,J),NR) .And. Corn(NR,J) /= B) then

						H = Corn(NR,J)

						if(.Not. IsInTheList(Points,PC,H)) then

							PC = PC + 1
							Points(PC) = H

						endif

						exit

					endif

				end do

				if(PC == 6) then !--------------- Four Quads MUST have formed a HEXAGON ------------
!Part 3:
					!================== Finding Quad Neibours surrounding hexagon ==================

					!--------------------------------- Specifying Ncd ------------------------------

                    Ncd = getNeibour(Dim,Corn,Neib,C,D,ME)

					!--------------------------------- Specifying Nde ------------------------------

                    Nde = getNeibour(Dim,Corn,Neib,E,D,NL)

					!--------------------------------- Specifying Nef ------------------------------

                    Nef = getNeibour(Dim,Corn,Neib,E,F,NL)

					!--------------------------------- Specifying Nfg ------------------------------

                    Nfg = getNeibour(Dim,Corn,Neib,G,F,NE)

					!--------------------------------- Specifying Ngh ------------------------------

                    Ngh = getNeibour(Dim,Corn,Neib,G,H,NR)

					!--------------------------------- Specifying Nhc ------------------------------

                    Nhc = getNeibour(Dim,Corn,Neib,C,H,NR)

					!======================== Calculating Nodes Valances =========================
!Part 4:
                    if(Ncd/=0 .And. Nde/=0 .And. Nef/=0 .And. Nfg/=0 .And. Ngh/=0 .And. Nhc/=0) then
                    
                        Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc)
                        Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)
                        Call GetSurroundingElements(Dim,Corn,Neib,E,NL,TElms,QElms,TEC,Ne1)
					    Call GetSurroundingElements(Dim,Corn,Neib,F,NE,TElms,QElms,TEC,Nf)
                        Call GetSurroundingElements(Dim,Corn,Neib,G,NE,TElms,QElms,TEC,Ng)
                        Call GetSurroundingElements(Dim,Corn,Neib,H,NR,TElms,QElms,TEC,Nh)

					    N1 = Nd + Ng
					    N2 = Nc + Nf
					
					    Q1 = (QuadQuality(Dim,X,Y,G,H,C,D) + QuadQuality(Dim,X,Y,E,F,G,D))/2
					    Q2 = (QuadQuality(Dim,X,Y,E,F,C,D) + QuadQuality(Dim,X,Y,C,F,G,H))/2

					    if(Ne1 >= 5 .And. Nh >= 5 .And. Nc ==4 .And. Nd == 4 .And. Nf == 4 .And. Ng == 4) then !----- Special Case -----
!Part 5:
						    if(SegmentJointIsPossible(Dim,X,Y,D,F,Points,PC) .And. SegmentJointIsPossible(Dim,X,Y,C,G,Points,PC)) then

							    done = .True.

							    !-------------------------------- Moving D to midpoint of DF ----------------------------

							    x_value = (X(D) + X(F))/2
							    y_value = (Y(D) + Y(F))/2

							    X(D) = x_value
							    Y(D) = y_value

							    !-------------------------------- Moving C to midpoint of CG ----------------------------

							    x_value = (X(C) + X(G))/2
							    y_value = (Y(C) + Y(G))/2

							    X(C) = x_value
							    Y(C) = y_value

							    !------------------------------------- Mapping F to D -----------------------------------

							    Call GetSurroundingElements(Dim,Corn,Neib,F,NL,TElms,QElms,TEC,QEC)
														
							    do J=1,QEC
							
								    Elm = QElms(J)

								    do K=1,4

									    if(Corn(Elm,K) == F) then

										    Corn(Elm,K) = D
										    exit

									    endif

								    end do

							    end do

							    !------------------------------------- Mapping G to C -----------------------------------

							    Call GetSurroundingElements(Dim,Corn,Neib,G,NR,TElms,QElms,TEC,QEC)
														
							    do J=1,QEC
							
								    Elm = QElms(J)

								    do K=1,4

									    if(Corn(Elm,K) == G) then

										    Corn(Elm,K) = C
										    exit

									    endif

								    end do

							    end do

							    !------------- Modifying Neibours of Neibours ---------

							    Call setNeibour(Dim,Corn,Neib,Nhc,H,C,Ngh)
							    Call setNeibour(Dim,Corn,Neib,Ngh,H,C,Nhc)

							    Call setNeibour(Dim,Corn,Neib,Ncd,C,D,Nfg)
							    Call setNeibour(Dim,Corn,Neib,Nfg,C,D,Ncd)
							
							    Call setNeibour(Dim,Corn,Neib,Nde,D,E,Nef)
							    Call setNeibour(Dim,Corn,Neib,Nef,D,E,Nde) 

							    !------------------ Deleting Elements -----------------

							    do J=2,4
								
								    Corn(NL,J) = Corn(NL,1)
								    Corn(ME,J) = Corn(ME,1)
								    Corn(NE,J) = Corn(NE,1)
								    Corn(NR,J) = Corn(NR,1)

							    end do

							    !======================== Smoothing Neibour Elements =========================
							
							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Ncd,J),Ncd)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nde,J),Nde)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nef,J),Nef)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nfg,J),Nfg)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Ngh,J),Ngh)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nhc,J),Nhc)

							    end do

						    endif !---------------------------- End of Segment Joint Checking ---------------------------
  

					    else

						    if((N1 < N2) .Or. (N1 == N2 .And. Q1 > Q2)) then !------------ State1: DG is used as new Segment -----------
!Part 6:
							    if(SegmentJointIsPossible(Dim,X,Y,D,G,Points,PC)) then

								    done = .True.

								    !====================== ABCD converts to CDGH ==========================

								    Corn(ME,Ai1) = G !---------------- A is mapped to G --------------------
								    Corn(ME,Bi1) = H !---------------- B is mapped to H --------------------

								    !====================== AFGB converts to DEFG ==========================

								    Corn(NE,Ai2) = E !---------------- A is mapped to E --------------------
								    Corn(NE,Bi2) = D !---------------- B is mapped to D --------------------

								    !======================== Modifying ME neibours ========================
                                
                                    Call setNeibour(Dim,Corn,Neib,ME,D,G,NE)  !--------- Setting NE as neibour on DG side ---------
                                    Call setNeibour(Dim,Corn,Neib,ME,H,G,Ngh) !--------- Setting Ngh as neibour on GH side --------
                                    Call setNeibour(Dim,Corn,Neib,ME,H,C,Nhc) !--------- Setting Nhc as neibour on HC side --------

								    !======================== Modifying NE neibours ========================

                                    Call setNeibour(Dim,Corn,Neib,NE,D,G,ME) !--------- Setting ME as neibour on DG side ----------
                                    Call setNeibour(Dim,Corn,Neib,NE,D,E,Nde) !--------- Setting Nde as neibour on DE side -------- 
                                    Call setNeibour(Dim,Corn,Neib,NE,F,E,Nef) !--------- Setting Nef as neibour on EF side ---------
                                
								    !=================== Modifying neibours of neibours ====================

								    !------------------------ Updating Nde Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nde,D,E,NE)

								    !------------------------ Updating Nef Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nef,F,E,NE)

								    !------------------------ Updating Ngh Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Ngh,G,H,ME)

								    !------------------------ Updating Nhc Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nhc,C,H,ME)

								    !====================== Deleting NL and NR elements ====================

                                    Call DeleteCell(Dim,NL,Corn)
                                    Call DeleteCell(Dim,NR,Corn)

							    endif
						
						    elseif((N1 > N2) .Or. (N1 == N2 .And. Q1 < Q2)) then !---------- State2: CF is used as new Segment -------------
!Part 7:
							    if(SegmentJointIsPossible(Dim,X,Y,C,F,Points,PC)) then

								    done = .True.

								    !====================== ABCD converts to EFCD ==========================

								    Corn(ME,Ai1) = E !---------------- A is mapped to E --------------------
								    Corn(ME,Bi1) = F !---------------- B is mapped to F --------------------

								    !====================== AFGB converts to CFGH ==========================

								    Corn(NE,Ai2) = C !---------------- A is mapped to E --------------------
								    Corn(NE,Bi2) = H !---------------- B is mapped to D --------------------

								    !======================== Modifying ME neibours ========================

                                    Call setNeibour(Dim,Corn,Neib,ME,C,F,NE) !--------- Setting NE as neibour on CF side ----------
                                    Call setNeibour(Dim,Corn,Neib,ME,E,F,Nef) !-------- Setting Nef as neibour on EF side ---------
                                    Call setNeibour(Dim,Corn,Neib,ME,E,D,Nde) !--------- Setting Nde as neibour on DE side --------

								    !======================== Modifying NE neibours ========================

                                    Call setNeibour(Dim,Corn,Neib,NE,C,F,ME) !--------- Setting ME as neibour on CF side ----------
                                    Call setNeibour(Dim,Corn,Neib,NE,C,H,Nhc) !--------- Setting Nhc as neibour on CH side --------
                                    Call setNeibour(Dim,Corn,Neib,NE,G,H,Ngh) !---------- Setting Ngh as neibour on GH side -------

								    !=================== Modifying neibours of neibours ====================

								    !------------------------ Updating Nde Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nde,D,E,ME)

								    !------------------------ Updating Nef Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nef,F,E,ME)

								    !------------------------ Updating Ngh Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Ngh,G,H,NE)

								    !------------------------ Updating Nhc Neibours ------------------------
								
                                    Call setNeibour(Dim,Corn,Neib,Nhc,C,H,NE)

								    !====================== Deleting NL and NR elements ====================

                                    Call DeleteCell(Dim,NL,Corn)
                                    Call DeleteCell(Dim,NR,Corn)

							    endif

						    endif	
!Part 8:
						    if(done) then
				
							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ME,J),ME)

							    end do

							    do J=1,4

								    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NE,J),NE)

							    end do

						    endif

                        endif
                    
                    endif

				endif !------------------------ End of Hexagon checking ----------------------------

			endif !---------------------- End of Segment Found condition --------------------------- 

		endif !-------------------------- End of Initial Condition -----------------------------

	endif !----------- End of Checking A not to be on the exterior boundary ----------------

	if(done) exit

end do
!===========================================================================================
End Subroutine SegmentElimination
!*********************************************************************************************
