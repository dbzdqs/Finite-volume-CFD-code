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
Subroutine OBS(Dim,NP,Corn,Neib,X,Y,Elm,V,DELTA,COINCIDENT_TOLERANCE,MOVE_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,NP,Corn,Neib,Elm,V,DELTA,COINCIDENT_TOLERANCE,MOVE_TOLERANCE
Intent(InOut)::X,Y

Integer,Parameter::Gx = 1
Integer,Parameter::Gy = 2
Real(8),Parameter::OBS_TOLERANCE = 0.1
Real(8),Parameter::TOL = 0.00001

Integer::Dim,NP,V,Elm,TEC,QEC,I,J,Iter,min_index,GC,State,E_min
Integer,Dimension(1:1000)::QElms,TElms
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::Node_Moved
Real(8)::DELTA,COINCIDENT_TOLERANCE,GAMMA,MOVE_TOLERANCE,QuadDistortionMetric,x_old,y_old,Mu_min,Mu,Mu_plus,Mu_min_new,Gdir_x,Gdir_y,G_norm,temp1,temp2,GAMMA_FACTOR
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:100)::PreDistortionMetrics,PostDistortionMetrics,GAMMA_Values
Real(8),Dimension(1:100,1:2)::Gradient
!===========================================================================================

!--->>> Reference: Canann, S. A., J. R. Tristano and M. L. Staten (1998), "An approach to combined Laplacian and
!--->>> optimization-based smoothing for triangular, quadrilateral and tetrahedral meshes"

GAMMA_FACTOR = COINCIDENT_TOLERANCE/10 

!Part 1:

Call GetSurroundingElements(Dim,Corn,Neib,V,Elm,TElms,QElms,TEC,QEC)

!--------------------------- Calculate Initial Distortion Metrics --------------------------

do I=1,QEC

	PreDistortionMetrics(I) = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE) 

end do

!Part 2:

!------------------------- Finding Minimum Distortion Metric(Mu_min) ----------------------- 

Mu_min = PreDistortionMetrics(1)
E_min = QElms(1)
min_index = 1

do I=2,QEC

	if(Mu_min > PreDistortionMetrics(I)) then
	
		Mu_min = PreDistortionMetrics(I)
		E_min = QElms(I)
		min_index = I

	endif

end do 

if(Mu_min <= OBS_TOLERANCE) then !-------------- Condition On Invoking OBS ------------------

	x_old = X(V)
	y_old = Y(V)
!Part 3:
	do Iter=1,2 !------------ Only Two Iterations are Sufficient ----------------------------

		!----------------------- Finding Gradient Vector per Element -------------------

		Call CalcGradientDirction(Dim,Corn,X,Y,V,DELTA,Gradient,PreDistortionMetrics,Mu_min,E_min,QElms,QEC,COINCIDENT_TOLERANCE)
		
		!------------------------ Finding Gradient Direction ---------------------------
		
		Gdir_x = Gradient(min_index,Gx)
		Gdir_y = Gradient(min_index,Gy)
		
		!---------------------------- Measuring Gdir norm ------------------------------

		G_norm = DSQRT(Gdir_x*Gdir_x + Gdir_y*Gdir_y)
		
		if(G_norm < 0.00001) then !---------- Gdir Must Not be too Small ---------------
		
			!------------------------- Finding Mu_min_plus -----------------------------
								
			Call getMU_min_plus(PreDistortionMetrics,QEC,Gradient,Gdir_x,Gdir_y)
		
        endif
            
        if(Gdir_x /= 0 .And. Gdir_y /= 0) then
            
		    G_norm = DSQRT(Gdir_x*Gdir_x + Gdir_y*Gdir_y)

		    !print *,'G: ',G_norm
		
		    !------------------------------- Calculating GAMMA ----------------------------- 

		    GC = 0

		    do I=1,QEC
		
			    temp1 = Gdir_x*Gdir_x + Gdir_y*Gdir_y
			    temp2 = Gdir_x*Gradient(I,Gx) + Gdir_y*Gradient(I,Gy)  

			    if(temp2 < 0) then
			
				    GC = GC + 1
			
				    Mu = PreDistortionMetrics(I)
				
				    GAMMA_Values(GC) = (Mu - Mu_min)/(temp1 - temp2) 

			    endif  

		    end do
		
		    if(GC  /= 0) then

			    GAMMA = GAMMA_Values(1) 

			    do I=1,GC
			
				    if(GAMMA > GAMMA_Values(I)) then

					    GAMMA = GAMMA_Values(I)	

				    endif

			    end do
		
            else

                Call CalcGamma(Dim,Corn,Neib,X,Y,V,Elm,Gdir_x,Gdir_y,GAMMA)	
		
		    endif
		
		    Node_Moved = .False.

		    !---------------------- Applying new Coordinate to V ---------------------------
		
		    X(V) = X(V) + GAMMA*Gdir_x
		    Y(V) = Y(V) + GAMMA*Gdir_y
!Part 4:
		    do I=1,4

			    !----------------------- Calculating new Distortion Metrics ----------------

			    do J=1,QEC

				    PostDistortionMetrics(J) = QuadDistortionMetric(Dim,Corn,X,Y,QElms(J),COINCIDENT_TOLERANCE) 		

			    end do

			    Mu_min_new = PostDistortionMetrics(1) 

			    do J=1,QEC
			
				    if(Mu_min_new > PostDistortionMetrics(J)) then

					    Mu_min_new = PostDistortionMetrics(J)	

				    endif

			    end do

			    if(Mu_min_new >= Mu_min + TOL) then

				    x_old = X(V)
				    y_old = Y(V)

				    Node_Moved = .True.

				    Call UpdateMetrics(Dim,Corn,X,Y,QElms,Mu_min,E_min,min_index,QEC,PreDistortionMetrics,COINCIDENT_TOLERANCE)

				    exit

			    else

				    GAMMA = GAMMA/2

				    X(V) = X(V) + GAMMA*Gdir_x
				    Y(V) = Y(V) + GAMMA*Gdir_y

			    endif

		    end do
		
		    if(.Not. Node_Moved) then
		
			    X(V) = x_old
			    Y(V) = y_old
                exit
		
            endif
            
        else
            
            exit
        
        endif

	end do !------------------------- End of 'Iterations' Loop -----------------------------

endif !-------------------------- End of Condition On Invoking OBS -------------------------

!===========================================================================================
End Subroutine OBS
!*********************************************************************************************
