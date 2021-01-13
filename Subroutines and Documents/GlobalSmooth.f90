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
Subroutine GlobalSmooth(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,NP,Corn,Neib,NBE,BFP
Intent(InOut)::X,Y

Integer,Parameter::NODE = 1
Integer,Parameter::ELEMENT = 2

Integer::Dim,NC,NP,NBE,Niter,E,I,V,MAXIMUM_Distance,PC,Num_Moved_Nodes
Integer,Dimension(1:Dim,1:2)::BFP,SmoothNode
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsOnExteriorBoundary,Inversion,QuadIsInverted
Logical,Dimension(1:Dim)::MovedByOBS,Deactive
Real(8)::DELTA,COINCIDENT_TOLERANCE,MOVE_TOLERANCE,x_old,y_old,Move,MoveDistance,LARGEST_DISTANCE_MOVE
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., J. R. Tristano and M. L. Staten (1998), "An approach to combined Laplacian and
!--->>> optimization-based smoothing for triangular, quadrilateral and tetrahedral meshes"

print *, '<<<<<<<<<<<<<<<<------ GLOBAL SMOOTH ------>>>>>>>>>>>>>>>'

!Part 1:

PC = 0 !--------------------- Number of Nodes to be processed ------------------------------

Call GetSmoothableNodes(Dim,NC,Corn,SmoothNode,PC)
Call InitLists(Dim,Deactive,MovedByOBS,PC)
Call CalculateConstants(Dim,NBE,BFP,X,Y,DELTA,COINCIDENT_TOLERANCE)

MOVE_TOLERANCE = COINCIDENT_TOLERANCE*0.01 

Niter = 1

DO  !--------------------------- Main Loop ------------------------------------

	Num_Moved_Nodes = 0
	LARGEST_DISTANCE_MOVE = 0

	DO I=1,PC !-------------------- Iterate over all smoothable Nodes ----------------------

		V = SmoothNode(I,NODE)	  !--------------- Selected Node to Process ----------------
		E = SmoothNode(I,ELEMENT) !-------------- Element of Selected Node -----------------

		IF(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,V)) THEN 

			x_old = X(V)
			y_old = Y(V)

			IF(.Not. MovedByOBS(V) .AND. .Not. Deactive(V)) THEN !----------- If V wasn't moved by OBS -------------
!Part 2:
				Call CLS(Dim,Corn,Neib,X,Y,E,V,COINCIDENT_TOLERANCE) !---- Constrained Laplacian Smoothing ----

				Move = MoveDistance(Dim,X,Y,V,x_old,y_old) 

				IF(Move < MOVE_TOLERANCE) THEN
				
					X(V) = x_old
					Y(V) = y_old

					Deactive(V) = .True.
					
				ELSE

					IF(LARGEST_DISTANCE_MOVE < Move) THEN

						LARGEST_DISTANCE_MOVE = Move 

					ENDIF
					
					Num_Moved_Nodes = Num_Moved_Nodes + 1 
					 
				ENDIF !---------------- End of Checking Move Distance ----------------------

			ENDIF !-------------- End of condition of If V movement by OBS -----------------

			IF(Niter >= 2) THEN
!Part 3:   
				!-> Invoke Optimization-based smoothing
                
				Call OBS(Dim,NP,Corn,Neib,X,Y,E,V,DELTA,COINCIDENT_TOLERANCE,MOVE_TOLERANCE)
    
				Move = MoveDistance(Dim,X,Y,V,x_old,y_old)
    
				IF(Move /= 0) THEN
					
					MovedByOBS(V) = .True.
    
					IF(LARGEST_DISTANCE_MOVE < Move) THEN
    
						LARGEST_DISTANCE_MOVE = Move 
    
					ENDIF
					
					Num_Moved_Nodes = Num_Moved_Nodes + 1
    
				ENDIF
   
			ENDIF

		ENDIF !---------------- End of Checking Mobility of Selected Node ------------------

    END DO !--------------------- Iterate over all smoothable Nodes ------------------------

	Niter = Niter + 1

    print*,'Number of Iterations: ',Niter,' : ',(LARGEST_DISTANCE_MOVE - (1.75*MOVE_TOLERANCE))

	IF( Num_Moved_Nodes == 0 .Or. LARGEST_DISTANCE_MOVE < (1.75*MOVE_TOLERANCE)) EXIT

END DO !----------------------------- End of Main Loop -------------------------------------
!===========================================================================================
End Subroutine GlobalSmooth
!*********************************************************************************************
