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
Subroutine QuadrilateralFormation(Dim,Corn,Neib,FrontElement,newQuad,NC,possible)
Implicit None
!===========================================================================================
Intent(In)::Dim,FrontElement,newQuad
Intent(Inout)::Corn,Neib,NC,possible

Integer::Dim,NC,FrontElement,DT_Count,DPL_Count,I,J,N1,N2,N3,CommonCorners,Q_Index1,Q_Index2,Q_Index3,DCI,CNDC,CCI,temp
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::DT,DPL
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::isQuadEdge,isQuadCorner,isInDeletingPointsList,haveEdgeInCommon,possible
!===========================================================================================
!Part 1:
!-------------------------------------------------------------- Defining new Quad cell -------------------------------------------
NC = NC + 1
Corn(NC,1) = newQuad(1)
Corn(NC,2) = newQuad(2)
Corn(NC,3) = newQuad(3)
Corn(NC,4) = newQuad(4)
!----------------------------------------------------- Initialization of neibours of new Quad ------------------------------------
Neib(NC,1) = 0
Neib(NC,2) = 0
Neib(NC,3) = 0
Neib(NC,4) = 0

!Part 2:
!---------------------------------------------------------- Gathering Deleting Triangles -----------------------------------------
DT_Count = 1
DT(DT_Count) = FrontElement
possible = .True.
Call GetDeletingTriangles(Dim,FrontElement,Corn,Neib,DT,DT_Count,newQuad,possible)

!_______________________________________________________ Start of Deletion Process __________________________________________

if(possible) then
    
    DPL_Count = 0
    do I=1,DT_Count
!Part 3:        
	    haveEdgeInCommon = .False.
	    !------------------------- Counting Number of Common Corners between Deleting triangle and new Quad --------------------------
	    CommonCorners = 0
	    do J=1,3
		    if(isQuadCorner(Corn(DT(I),J),newQuad)) then
			    CommonCorners = CommonCorners + 1	
		    endif
	    end do
    !__________________________________________ Situation in which triangles have no common edge with Quad ___________________________
!Part 4:
	    if(CommonCorners == 0) then
	    !-------------------- Deleting triangle Cell and Adding its corners to DPL (Deleting Points List) ----------------------------	
		    do J=1,3
			    if(.Not. isInDeletingPointsList(DPL,DPL_Count,Corn(DT(I),J))) then
				    DPL_Count = DPL_Count + 1
				    DPL(DPL_Count) = Corn(DT(I),J)
			    endif	
		    end do
    !______________________________________ Situation in which triangles have one or two common edge with Quad ______________________________________________________

	    elseif(CommonCorners == 2 .Or. CommonCorners == 3) then
		
		    !-------------------- Situation in which current triangle and the new Quad have two edges in common ----------------------
		
		    if(CommonCorners == 3) then
!Part 5:			
			    !-------------------- Specifying indices of corners of new Quad that are common with deleting triangle ---------------  
			    do J=1,4
				    if(newQuad(J) == Corn(DT(I),1)) Q_Index1 = J
				    if(newQuad(J) == Corn(DT(I),2)) Q_Index2 = J
				    if(newQuad(J) == Corn(DT(I),3)) Q_Index3 = J
			    end do

			    if(isQuadEdge(Corn(DT(I),1),Corn(DT(I),2),newQuad) .And. isQuadEdge(Corn(DT(I),2),Corn(DT(I),3),newQuad)) then

				    if(Neib(DT(I),1) /= 0) then !------------------------- Processing Edge 1->2 --------------------------------------
					    !------------------ Specifying neibour of neibour of deleting triangle as neibour of new Quad ----------------
					    N1 = Neib(DT(I),1)
					    if(Q_Index2 < Q_Index3) then
						    if(Q_Index3 - Q_Index2 == 1) then
							    Neib(NC,Q_Index2) = N1
						    else
							    Neib(NC,Q_Index3) = N1
						    endif
					    elseif(Q_Index2 > Q_Index3) then
						    if(Q_Index2 - Q_Index3 == 1) then
							    Neib(NC,Q_Index3) = N1
						    else
							    Neib(NC,Q_Index2) = N1
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N1,J) == DT(I)) then
							    Neib(N1,J) = NC 
						    endif
					    end do
				    endif

				    if(Neib(DT(I),3) /= 0) then !------------------------- Processing Edge 2->3 --------------------------------------
				    !-------------------- Specifying neibour of neibour of deleting triangle as neibour of new Quad ------------------
					    N3 = Neib(DT(I),3)
					    if(Q_Index1 < Q_Index2) then
						    if(Q_Index2 - Q_Index1 == 1) then
							    Neib(NC,Q_Index1) = N3
						    else
							    Neib(NC,Q_Index2) = N3
						    endif
					    elseif(Q_Index1 > Q_Index2) then
						    if(Q_Index1 - Q_Index2 == 1) then
							    Neib(NC,Q_Index2) = N3
						    else
							    Neib(NC,Q_Index1) = N3
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N3,J) == DT(I)) then
							    Neib(N3,J) = NC
						    endif
					    end do
				    endif

			    elseif(isQuadEdge(Corn(DT(I),2),Corn(DT(I),3),newQuad) .And. isQuadEdge(Corn(DT(I),3),Corn(DT(I),1),newQuad)) then

				    if(Neib(DT(I),1) /= 0) then !------------------------- Processing Edge 3->1 --------------------------------------
					    !------------------ Specifying neibour of neibour of deleting triangle as neibour of new Quad ----------------
					    N1 = Neib(DT(I),1)
					    if(Q_Index2 < Q_Index3) then
						    if(Q_Index3 - Q_Index2 == 1) then
							    Neib(NC,Q_Index2) = N1
						    else
							    Neib(NC,Q_Index3) = N1
						    endif
					    elseif(Q_Index2 > Q_Index3) then
						    if(Q_Index2 - Q_Index3 == 1) then
							    Neib(NC,Q_Index3) = N1
						    else
							    Neib(NC,Q_Index2) = N1
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N1,J) == DT(I)) then
							    Neib(N1,J) = NC 
						    endif
					    end do
				    endif

				    if(Neib(DT(I),2) /= 0) then !------------------------- Processing Edge 2->3 --------------------------------------
				    !-------------------- Specifying neibour of neibour of deleting triangle as neibour of new Quad ------------------
					    N2 = Neib(DT(I),2)
					    if(Q_Index1 < Q_Index3) then
						    if(Q_Index3 - Q_Index1 == 1) then
							    Neib(NC,Q_Index1) = N2
						    else
							    Neib(NC,Q_Index3) = N2
						    endif
					    elseif(Q_Index1 > Q_Index3) then
						    if(Q_Index1 - Q_Index3 == 1) then
							    Neib(NC,Q_Index3) = N2
						    else
							    Neib(NC,Q_Index1) = N2
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N2,J) == DT(I)) then
							    Neib(N2,J) = NC
						    endif
					    end do
				    endif

			    elseif(isQuadEdge(Corn(DT(I),1),Corn(DT(I),2),newQuad) .And. isQuadEdge(Corn(DT(I),3),Corn(DT(I),1),newQuad)) then
																				 
				    if(Neib(DT(I),2) /= 0) then !------------------------- Processing Edge 1->2 --------------------------------------
					    !------------------ Specifying neibour of neibour of deleting triangle as neibour of new Quad ----------------
					    N2 = Neib(DT(I),2)
					    if(Q_Index1 < Q_Index3) then
						    if(Q_Index3 - Q_Index1 == 1) then
							    Neib(NC,Q_Index1) = N2
						    else
							    Neib(NC,Q_Index3) = N2
						    endif
					    elseif(Q_Index1 > Q_Index3) then
						    if(Q_Index1 - Q_Index3 == 1) then
							    Neib(NC,Q_Index3) = N2
						    else
							    Neib(NC,Q_Index1) = N2
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N2,J) == DT(I)) then
							    Neib(N2,J) = NC
						    endif
					    end do
				    endif

				    if(Neib(DT(I),3) /= 0) then !------------------------- Processing Edge 3->1 --------------------------------------
				    !-------------------- Specifying neibour of neibour of deleting triangle as neibour of new Quad ------------------
					    N3 = Neib(DT(I),3)
					    if(Q_Index1 < Q_Index2) then
						    if(Q_Index2 - Q_Index1 == 1) then
							    Neib(NC,Q_Index1) = N3
						    else
							    Neib(NC,Q_Index2) = N3
						    endif
					    elseif(Q_Index1 > Q_Index2) then
						    if(Q_Index1 - Q_Index2 == 1) then
							    Neib(NC,Q_Index2) = N3
						    else
							    Neib(NC,Q_Index1) = N3
						    endif
					    endif
					    !------------------------- Specifying new Quad as neibour of neibour of deleting triangle --------------------
					    do J=1,4
						    if(Neib(N3,J) == DT(I)) then
							    Neib(N3,J) = NC
						    endif
					    end do
				    endif
			    endif
		    !-------------------- Situation in which current triangle and the new Quad have one edge in common -----------------------
            elseif(CommonCorners == 2) then
!Part 6:                
			    !--------------------------- Determining Index of Deleting Corner of current triangle --------------------------------
			    do J=1,3
				    if(.Not. (isQuadCorner(Corn(DT(I),J),newQuad))) then
					    DCI = J !--------------------------------------- DCI (Deleting Corner Index) ---------------------------------
				    endif
			    end do
			    !------------------------------------- Checking if triangle and Quad have an Edge in common --------------------------
			    if(isQuadEdge(Corn(DT(I),1),Corn(DT(I),2),newQuad)) then
				    haveEdgeInCommon = .True.
			    elseif(isQuadEdge(Corn(DT(I),2),Corn(DT(I),3),newQuad)) then
				    haveEdgeInCommon = .True.
			    elseif(isQuadEdge(Corn(DT(I),3),Corn(DT(I),1),newQuad)) then
				    haveEdgeInCommon = .True.
			    endif

			    if(haveEdgeInCommon) then

				    !-------------------------------- Specifying new Quad as neibour of neibour of deleting triangle -----------------
				    if(Neib(DT(I),DCI) /= 0) then
					    CNDC = Neib(DT(I),DCI) !----------------- CNDC (Corresponding Neibour of Deleting Corner) -------------------- 
					    do J=1,4
						    if(Neib(CNDC,J) == DT(I)) then
							    Neib(CNDC,J) = NC
						    endif
					    end do
				    endif
				    !-------------------- Specifying neibour of neibour of deleting triangle as neibour of new Quad ------------------
				    if(DCI == 1) then
					    !---------------- Specifying indices of corners of new Quad that are common with deleting triangle -----------
					    do J=1,4
						    if(newQuad(J) == Corn(DT(I),2)) Q_Index2 = J
						    if(newQuad(J) == Corn(DT(I),3)) Q_Index3 = J
					    end do
					    !------------------------------------------- Defining Neibour of new Quad ------------------------------------
					    if(Neib(DT(I),1) /= 0) then
						    N1 = Neib(DT(I),1)
						    if(Q_Index2 < Q_Index3) then
							    if(Q_Index3 - Q_Index2 == 1) then 
								    Neib(NC,Q_Index2) = N1
							    else
								    Neib(NC,Q_Index3) = N1
							    endif
						    elseif(Q_Index2 > Q_Index3) then
							    if(Q_Index2 - Q_Index3 == 1) then
								    Neib(NC,Q_Index3) = N1
							    else
								    Neib(NC,Q_Index2) = N1
							    endif
						    endif
					    endif
				    elseif(DCI == 2) then
					    !---------------- Specifying indices of corners of new Quad that are common with deleting triangle -----------
					    do J=1,4
						    if(newQuad(J) == Corn(DT(I),1)) Q_Index1 = J
						    if(newQuad(J) == Corn(DT(I),3)) Q_Index3 = J
					    end do
					    !------------------------------------------- Defining Neibour of new Quad ------------------------------------
					    if(Neib(DT(I),2) /= 0) then
						    N2 = Neib(DT(I),2)
						    if(Q_Index1 < Q_Index3) then
							    temp = Q_Index3 - Q_Index1
							    if(temp == 1) then
								    Neib(NC,Q_Index1) = N2
							    else
								    Neib(NC,Q_Index3) = N2
							    endif
						    elseif(Q_Index1 > Q_Index3) then
							    temp = Q_Index1 - Q_Index3 
							    if(temp == 1) then
								    Neib(NC,Q_Index3) = N2
							    else
								    Neib(NC,Q_Index1) = N2
							    endif
						    endif
					    endif
				    elseif(DCI == 3) then
					    !---------------- Specifying indices of corners of new Quad that are common with deleting triangle -----------
					    do J=1,4
						    if(newQuad(J) == Corn(DT(I),1)) Q_Index1 = J
						    if(newQuad(J) == Corn(DT(I),2)) Q_Index2 = J
					    end do
					    !------------------------------------------- Defining Neibour of new Quad ------------------------------------
					    if(Neib(DT(I),3) /= 0) then
						    N3 = Neib(DT(I),3)
						    if(Q_Index1 < Q_Index2) then
							    if(Q_Index2 - Q_Index1 == 1) then
								    Neib(NC,Q_Index1) = N3
							    else 
								    Neib(NC,Q_Index2) = N3
							    endif
						    elseif(Q_Index1 > Q_Index2) then
							    if(Q_Index1 - Q_Index2 == 1) then
								    Neib(NC,Q_Index2) = N3
							    else
								    Neib(NC,Q_Index1) = N3
							    endif
						    endif
					    endif
				    endif

			    endif
			    !--------------------------- Adding triangle's extra points to DPL (Deleting Points List) ----------------------------
			    do J=1,3
				    if(J == DCI) then
					    if(.Not. isInDeletingPointsList(DPL,DPL_Count,Corn(DT(I),J))) then
						    DPL_Count = DPL_Count + 1
						    DPL(DPL_Count) = Corn(DT(I),J)
					    endif		
				    endif
			    end do
				
		    endif
	    !------------------------- Situation in which Deleting triangle and new Quad have one point in common ------------------------
        elseif(CommonCorners == 1) then
!Part 7:            
		    !---------------------------------------- Finding Index of the Common Corner ---------------------------------------------
		    do J=1,3
			    if(isQuadCorner(Corn(DT(I),J),newQuad)) then
				    CCI = J !------------------------------------ CCI (Common Corner Index) ------------------------------------------
			    endif
		    end do
		    !------------------- Deleting triangle cell and adding its extra points to DPL (Deleting Points List) --------------------
		    do J=1,3
			    if(J /= CCI) then
				    if(.Not. isInDeletingPointsList(DPL,DPL_Count,Corn(DT(I),J))) then
					    DPL_Count = DPL_Count + 1
					    DPL(DPL_Count) = Corn(DT(I),J)
				    endif	
			    endif
		    end do			
	    endif
    end do
!Part 8:    
    !-------------------------------------------- Deleting Extra Cells and Their Neibours --------------------------------------------
    do I=1,DT_Count
	    print *,'DT',I,':',DT(I)
	    Call DeleteCell(Dim,DT(I),Corn)
    end do
else
    Call DeleteCell(Dim,NC,Corn)    
endif
!===========================================================================================
End Subroutine QuadrilateralFormation
!*********************************************************************************************
