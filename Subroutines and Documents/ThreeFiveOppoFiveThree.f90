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
Subroutine ThreeFiveOppoFiveThree(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,ME,done) !--->>> Pattern in Fig 32 of reference
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,TEC,ME,newElement,newN,I,J,A,Ai,B,Bi,C,Ci,D,Di,E,F,Na,Nb,Nc1,Nd,Ne,Nf,Nab,Nbc,Ncd,Nda,Elm,getNeibour,Cnt
Integer,Dimension(1:1000)::QElms,TElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,done1,areAdjacent,IsDiagonalInsideQuad,IsOnExteriorBoundary
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do I=1,4

	done = .False.
	done1 = .False.

	Ai = I
	A = Corn(ME,Ai)
!Part 1:    
    if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,A)) then !--- checking A not to be on exterior boundary ---
	
	    Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElms,TEC,Na)

	    if(Na == 5) then !------------------------ Condition on A Valance ----------------------

		    Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

            if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,C)) then 
            
		        Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc1)

		        if(Nc1 == 5) then !--------------------- Condition on C Valance --------------------

			        Call getNextCorner(Dim,Corn,ME,A,B,Bi)
		
			        E = B
			        Elm = getNeibour(Dim,Corn,Neib,A,B,ME)
			        Cnt = 0

			        do

				        do J=1,4

					        if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= E) then

						        Cnt = Cnt + 1
							
						        E = Corn(Elm,J)

						        if(Cnt /= 2) then

							        Elm = getNeibour(Dim,Corn,Neib,A,E,Elm)

						        endif

						        exit

					        endif

				        end do

				        if(Cnt == 2) exit

                    end do

                    if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,E)) then
                    
			            Call GetSurroundingElements(Dim,Corn,Neib,E,Elm,TElms,QElms,TEC,Ne)

			            if(Ne == 3) then !------------------- Condition on E Valance -------------------

				            F = B
				            Elm = getNeibour(Dim,Corn,Neib,C,B,ME)
				            Cnt = 0

				            do

					            do J=1,4

						            if(areAdjacent(Dim,Corn,C,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= F) then

							            Cnt = Cnt + 1
								
							            F = Corn(Elm,J)

							            if(Cnt /= 2) then

								            Elm = getNeibour(Dim,Corn,Neib,C,F,Elm)

							            endif

							            exit

						            endif

					            end do

					            if(Cnt == 2) exit

                            end do
				
                            if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,F)) then
                            
				                Call GetSurroundingElements(Dim,Corn,Neib,F,Elm,TElms,QElms,TEC,Nf)
				
				                if(Nf == 3) then !----------------- Condition on F Valance -----------------
!Part 2:				
					                if(IsDiagonalInsideQuad(Dim,Corn,X,Y,ME,A,C)) then

						                Call getPrevCorner(Dim,Corn,ME,A,D,Di)

						                !------------------------- Specifying ME Neibours ----------------------

						                Nab = getNeibour(Dim,Corn,Neib,A,B,ME)
						                Nbc = getNeibour(Dim,Corn,Neib,B,C,ME)
						                Ncd = getNeibour(Dim,Corn,Neib,C,D,ME)
						                Nda = getNeibour(Dim,Corn,Neib,D,A,ME)
!Part 3:
						                !---------------------- Adding a Two-Edged Node to ME ------------------
						
						                NP = NP + 1
						                newN = NP
						
						                X(newN) = (X(A)+X(C))/2
						                Y(newN) = (Y(A)+Y(C))/2
						
						                !---------------------------- Modifying ME -----------------------------
						
						                Corn(ME,Di) = newN

						                !-------------------- Adding new Element (based on ME) -----------------

						                NC = NC + 1
						                newElement = NC

						                Corn(newElement,Ai) = A
						                Corn(newElement,Bi) = newN
						                Corn(newElement,Ci) = C
						                Corn(newElement,Di) = D

						                !-------------------------- Modifying ME Neibours ----------------------

						                Call setNeibour(Dim,Corn,Neib,ME,A,newN,newElement)
						                Call setNeibour(Dim,Corn,Neib,ME,C,newN,newElement)

						                !-------------------- Specifying newElement Neibours -------------------

						                Call setNeibour(Dim,Corn,Neib,newElement,A,newN,ME)
						                Call setNeibour(Dim,Corn,Neib,newElement,C,newN,ME)
						                Call setNeibour(Dim,Corn,Neib,newElement,C,D,Ncd)
						                Call setNeibour(Dim,Corn,Neib,newElement,D,A,Nda)

						                !--------------------- Modifying Neibours of Neibours ------------------

						                Call setNeibour(Dim,Corn,Neib,Ncd,C,D,newElement)
						                Call setNeibour(Dim,Corn,Neib,Nda,D,A,newElement)

						                !-------------------------- Resolving Pattern ---------------------------
!Part 4:
						                Call ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,newN,A,E,ME,done)
						                Call ElementOpenOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,newN,C,F,ME,done1)

						                if(done .And. done1) then

							                exit

						                else
!Part 5:
							                done = .False.

							                Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,ME,done1)	

						                endif
					
					                endif 

                                endif !------------------- End of Condition on F Valance -------------------
                            
                            endif !------------------ End of checking F on exterior bounadary ----------

                        endif !--------------------- End of Condition on E Valance ---------------------
                    
                    endif !--------------------- End of Checking E on exterior boundary ------------

                endif !----------------------- End of Condition on C Valance -----------------------
            
            endif !-------------- End of checking C not to be on exterior boundary -------------

        endif !-------------------------- End of Condition on A Valance ------------------------
    
    endif !------------------------ End of checking A not to be on exterior boundary -------

end do !------------------------- End of Main Loop (I Loop) --------------------------------
!===========================================================================================
End Subroutine ThreeFiveOppoFiveThree 
!*********************************************************************************************
