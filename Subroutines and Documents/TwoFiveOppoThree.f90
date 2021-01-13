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
Subroutine TwoFiveOppoThree(Dim,Corn,Neib,X,Y,NBE,BFP,ME,done) !--->>> Pattern in Fig 31 of reference
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::Corn,Neib,X,Y,done

Integer::Dim,NBE,TEC,ME,NE1,NE2,NL,NR,K,L,A,Ai,B,Bi,C,Ci,D,Di,E,F,G,H,I,J,Na,Nb,Nc,Nd,Ne,Nf,Ng,Nh,Ni,Nj,Pi,Nbg,Nde,getNeibour,PC
Integer,Dimension(1:1000)::QElms,TElms,PList
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,areAdjacent,SegmentJointIsPossible,QuadIsInverted
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do K=1,4

	done = .False.

	PC = 0 !------------- A Counter for number of boundary points of pattern ---------------

	Ai = K
	A = Corn(ME,Ai)
!Part 1:	
	Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElms,TEC,Na)

	if(Na == 3) then !------------------------ Condition on A Valance ----------------------

		Call getNextCorner(Dim,Corn,ME,A,B,Bi)
		Call getOppoCorner(Dim,Corn,ME,A,C,Ci)
		Call getPrevCorner(Dim,Corn,ME,A,D,Di)

		Call addToList(PList,PC,B)
		Call addToList(PList,PC,C)
		Call addToList(PList,PC,D) 

		Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElms,TEC,Nb)
		Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc)
		Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)

		if(Nb == 4 .And. Nc == 4 .And. Nd == 4) then !---- Condition of other ME corners ---

			NE1 = getNeibour(Dim,Corn,Neib,A,B,ME)
			NE2 = getNeibour(Dim,Corn,Neib,A,D,ME)
			
            if(NE1 /= 0 .And. NE2 /= 0) then !---------------------- Checking if NE1,NE2 exists ------------------
                
			    Call getOppoCorner(Dim,Corn,NE1,A,G,Pi)
			    Call getOppoCorner(Dim,Corn,NE2,A,E,Pi)

			    Call addToList(PList,PC,E)
			
			    Call GetSurroundingElements(Dim,Corn,Neib,G,NE1,TElms,QElms,TEC,Ng)
			    Call GetSurroundingElements(Dim,Corn,Neib,E,NE2,TElms,QElms,TEC,Ne)
			
			    if(Ng == 5 .And. Ne == 5) then

				    !----------- Finding 'F' (other common corner between NE1 and NE2) ---------

				    do L=1,4

					    if(areAdjacent(Dim,Corn,A,Corn(NE1,L),NE1) .And. Corn(NE1,L) /= B) then

						    F = Corn(NE1,L)
						    exit

					    endif

				    end do

				    Call GetSurroundingElements(Dim,Corn,Neib,F,NE1,TElms,QElms,TEC,Nf)

				    if(Nf == 4) then

					    NL = getNeibour(Dim,Corn,Neib,F,G,NE1)
					    NR = getNeibour(Dim,Corn,Neib,F,E,NE2)
					
					    Call getOppoCorner(Dim,Corn,NL,G,I,Pi)
					    Call getOppoCorner(Dim,Corn,NL,F,H,Pi)
					    Call getOppoCorner(Dim,Corn,NR,F,J,Pi)

					    Call addToList(PList,PC,J)
					    Call addToList(PList,PC,I)
					    Call addToList(PList,PC,H)
					    Call addToList(PList,PC,G)
					
					    Call GetSurroundingElements(Dim,Corn,Neib,H,NL,TElms,QElms,TEC,Nh)
					    Call GetSurroundingElements(Dim,Corn,Neib,I,NL,TElms,QElms,TEC,Ni)
					    Call GetSurroundingElements(Dim,Corn,Neib,J,NR,TElms,QElms,TEC,Nj)
					
					    if(Nh == 4 .And. Ni == 4 .And. Nj == 4) then !---- H,I,J conditions ----

						    if(SegmentJointIsPossible(Dim,X,Y,B,I,PList,PC) .And. SegmentJointIsPossible(Dim,X,Y,D,I,PList,PC)) then
!Part 2:
							    done = .True.

							    !---------------- Finding Neccessary Neibours -------------------

							    Nbg = getNeibour(Dim,Corn,Neib,B,G,NE1)
							    Nde = getNeibour(Dim,Corn,Neib,D,E,NE2) 

							    !================= Modifying ME,NL,NR Elements ==================

							    Corn(ME,Ai) = I !--------------- Mapping A to I -----------------

							    do L=1,4

								    !------------------- Mapping F to B in NL -------------------

								    if(Corn(NL,L) == F) then

									    Corn(NL,L) = B

								    endif

								    !------------------- Mapping F to D in NR -------------------

								    if(Corn(NR,L) == F) then

									    Corn(NR,L) = D

								    endif

							    end do

							    !------------------- Checking Inversion of ME,NL,NR ------------

							    if(QuadIsInverted(Dim,Corn,X,Y,ME) .Or. QuadIsInverted(Dim,Corn,X,Y,NL) .Or. QuadIsInverted(Dim,Corn,X,Y,NR)) then
!Part 3:
								    done = .False.

								    Corn(ME,Ai) = A !----------- Reversing A Mapping -----------

								    do L=1,4

									    !------------ Reversing Mapping F to B in NL -----------

									    if(Corn(NL,L) == B) then

										    Corn(NL,L) = F

									    endif

									    !----------- Reversing A Mapping F to D in NR ----------

									    if(Corn(NR,L) == D) then

										    Corn(NR,L) = F

									    endif

								    end do

							    else
!Part 4:
								    !==================== Modifying ME Neibours ====================

								    Call setNeibour(Dim,Corn,Neib,ME,I,B,NL)
								    Call setNeibour(Dim,Corn,Neib,ME,D,I,NR)

								    !==================== Modifying NL Neibours ====================

								    Call setNeibour(Dim,Corn,Neib,NL,I,B,ME)
								    Call setNeibour(Dim,Corn,Neib,NL,G,B,Nbg)

								    !==================== Modifying NR Neibours ====================

								    Call setNeibour(Dim,Corn,Neib,NR,I,D,ME)
								    Call setNeibour(Dim,Corn,Neib,NR,D,E,Nde)

								    !================ Modifying Neibours of Neibours ===============

								    Call setNeibour(Dim,Corn,Neib,Nbg,B,G,NL)
								    Call setNeibour(Dim,Corn,Neib,Nde,D,E,NR)

								    !================== Deleting NE1,NE2 Elements ==================

								    do L=2,4

									    Corn(NE1,L) = Corn(NE1,1)
									    Corn(NE2,L) = Corn(NE2,1)

								    end do

								    !------------------------ Smoothing Elements -----------------------

								    do L=1,4

									    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(ME,L),ME)

								    end do

								    do L=1,4

									    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NL,L),NL)

								    end do

								    do L=1,4

									    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NR,L),NR)

								    end do

								    exit
							
							    endif	   

						    endif !---------- End of Segment Joint Checkings -------------------
					
					    endif !-------------------- End of H,I,J conditions --------------------	

				    endif !----------------- End of Condition on F Valance ---------------------

                endif !---------------- End of conditions on G and E Valence ------------------- 

            endif !----------------- End of Checking if NE1,NE2 exists ---------------------
        
		endif !------------------- End of Condition of other ME corners --------------------
	
	endif !------------------------- End of Condition on A Valance ------------------------- 

end do !-------------------------- End of Main Loop ('K' Loop) -----------------------------
!===========================================================================================
End Subroutine TwoFiveOppoThree
!*********************************************************************************************
