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
Subroutine DefineSideEdges(Dim,NC,NP,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,Angles,X,Y,current_Level,CF,of_State) !Tested
Implicit None
!===========================================================================================
Intent(In)::Dim,NBC,BFP,current_Level,CF,of_State
Intent(Inout)::Corn,Neib,X,Y,States,Angles,Fronts,FrontEdges,NC,NP

Integer,Parameter::STATUS_ZERO = 0
Integer,Parameter::STATUS_ONE = 1
Integer,Parameter::STATUS_TWO = 2
Integer,Parameter::STATUS_THREE = 3

Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4
Integer,Parameter::PF=2	!---- PF(Processed Fronts)

Integer::Dim,Fronts,NP,NC,NBC,MFEC,NFEC,NFE,NFI,NF_Vertex_Left,NF_Vertex_Right,CF,NK,NV,ME,MQE,NQE,BQE,TEC,QEC,current_Level,J,of_State
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::isNotChevron,TopEdgeRecoveryIsPossible,SpecialCasePossible,QFPossible,isThreeEdgeLoop
Real(8)::ANR,ANL,pi,MFN,NFN,Ratio,GetNorm,alpha1,alpha2
Real(8),Dimension(1:Dim)::X,Y
Real(8),Dimension(1:Dim,1:2)::Angles
!===========================================================================================
pi = 3.14159265358979d0
alpha1 = 0.04*pi
alpha2 = 0.09*pi
!Part 1:
Call CheckThreeEdgeLoop(Dim,Corn,Neib,FrontEdges,States,Fronts,CF,isThreeEdgeLoop) 

if(.Not. isThreeEdgeLoop) then

    Call InitQuad(newQuad) !---- Buffer Initialization in order to Avoid First Chance Exception Occurance ----

    Select Case(of_State)
    ! _________________________________ Processing Fronts With the Highest(First) Priority (State 11) ________________________________
	
        Case(STATUS_THREE)
!Part 2:
            print *, 'ANL: ',Angles(CF,LeftVertex)
            print *, 'ANR: ',Angles(CF,RightVertex)
                        
            NK = FrontEdges(CF,LeftVertex)
            NV = FrontEdges(CF,RightVertex)
                        
            newQuad(3)=NK								            ! Saving point number of third corner of new Quad
                        
            Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,1,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Left)
                        
            NK = FrontEdges(CF,RightVertex)
            NV = FrontEdges(CF,LeftVertex)
                        
            newQuad(4)=NK								            ! Saving point number of fourth corner of new Quad
                        
            Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,2,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Right)
                        
            newQuad(1)=NF_Vertex_Right							 ! Saving point number of first corner of new Quad
            newQuad(2)=NF_Vertex_Left								 ! Saving point number of second corner of new Quad
						
            print *,'F:',CF,'&',':','newQuad3: ',newQuad(1),newQuad(2),newQuad(3),newQuad(4)
                        
            Call TopEdgeRecovery(Dim,NC,newQuad(1),newQuad(2),Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
						
            if(TopEdgeRecoveryIsPossible) then

	            Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
                            
                if(QFPossible) then
                                
		            Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
		            Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
		            Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                            
                else
                                
                    print *,'---------->>>>>>>> Quad Formation Impossible! <<<<<<<<----------'
                                
                    Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                                
                endif

            else
							
	            Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	            Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	            Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	            Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	            Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

            endif
		
    ! ___________________________________________ End of Processing Fronts With the Highest(First) Priority _______________________________

    ! ______________________________________________ Processing Fronts With Second Priority (State 01) _______________________________

	    Case(STATUS_ONE)

	    !-------------------- Situation in which one side edge on the left should be defined -------------------------
!Part 3:
		    !------------------------------ Initializing Base Edge of new Quadrilateral ----------------------------------

		    newQuad(3)=FrontEdges(CF,LeftVertex)                   ! Saving point number of third corner of new Quad
		    newQuad(4)=FrontEdges(CF,RightVertex)                  ! Saving point number of fourth corner of new Quad

		    !------------ In this loop we find the corner of the Main Front that is not part of its Front Edge ------------

		    do J=1,3									   
			    if(Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,LeftVertex) .And. Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,RightVertex)) then
				    MFEC=Corn(FrontEdges(CF,Element),J)          ! MFEC(Main Front Element Corner)
				    MQE = Neib(FrontEdges(CF,Element),J)		 ! MQE(Main Quad Element of Current Front)
			    endif								   
		    end do

		    NK = FrontEdges(CF,RightVertex)
		    NV = FrontEdges(CF,LeftVertex)				   ! The other corner of the current front edge called NV
		    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
		    ANR = Angles(CF,2)

		    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,2,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Right)

		    do J=1,3
			    if(Corn(NFE,J) == NFEC) then
				    NQE = Neib(NFE,J)
				    exit 
			    endif
		    end do

		    MFN = GetNorm(Dim,NK,NV,X,Y)                   ! Main Front Norm
		    NFN = GetNorm(Dim,NK,NF_Vertex_Right,X,Y) 	   ! Neibour Front Norm
		    Ratio = DMAX1(MFN,NFN)/DMIN1(MFN,NFN)

		    !----------------------------------- Finding Quad Element of Bigger Front -----------------------------------
		    if(MFN > NFN) then
			    BQE = MQE
		    else
			    BQE = NQE
		    endif

		    Call GetSurroundingElements(Dim,Corn,Neib,NK,ME,TElms,QElms,TEC,QEC)

		    if((QEC>=5 .And. ANR<Alpha1) .Or. (QEC<5 .And. ANR<Alpha2)) then !------------------------------ Seaming ----------------------------------
!Part 4:
			    if(Ratio <= 2.5) then
   
				    print *,'------->>>>>>>>>>>> Seam Operation In STATUS_ONE <<<<<<<<<<<<-------'
				
				    Call Seam(Dim,NC,NP,NK,NV,NF_Vertex_Right,CF,NFI,Corn,Neib,X,Y,FrontEdges,Fronts,States,current_Level,1,SpecialCasePossible)
				
				    if(SpecialCasePossible) then
   
					    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
   
				    else
                    
					    print *,'---------->>>>>>>>>> Seam Impossible!!! <<<<<<<<<<<<<-----------'
                    
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
   
				    endif

			    elseif(Ratio > 2.5) then

				    if(BQE/=0 .And. isNotChevron(Dim,Corn,X,Y,BQE)) then

					    print *,'------->>>>>>>>>>>> Transition Seam in STATUS_ONE <<<<<<<<<<<<-------'
!Part 5:                    
					    Call TransitionSeam(Dim,NC,NP,NBC,NK,NV,NF_Vertex_Right,CF,NFI,ME,NFE,MFEC,NFEC,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,SpecialCasePossible)
					
                        if(SpecialCasePossible) then
                        
                            Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                        
                        else
                    
                            print *,'------->>>>>>>>>>>> Transition Seam is Impossible <<<<<<<<<<<<-------'
                        
                            Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                        Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                        
                        endif

                    else
					
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

				    endif

			    endif

		    else

			    if(Ratio >= 2.5 .And. BQE/=0 .And. isNotChevron(Dim,Corn,X,Y,BQE)) then !--------------------- Transition Split ------------------------------

				    print *,'------->>>>>>>>>>>> Transition Split in STATUS_ONE <<<<<<<<<<<<-------'
!Part 6:                
				    Call TransitionSplit(Dim,NC,NP,NBC,NK,NV,NF_Vertex_Right,ME,NFE,MFEC,NFEC,CF,NFI,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,SpecialCasePossible)
                
                    if(SpecialCasePossible) then
                    
				        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                    else
                    
                        print *,'------->>>>>>>>>>>> Transition Split is Impossible <<<<<<<<<<<<-------'
                    
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                    endif

			    else
!Part 7:
				    newQuad(1)=NF_Vertex_Right							 ! Saving point number of first corner of new Quad

				    !----------------------------------- Defining Side Edge on Left Vertex of Front Edge------------------------------------------
			
				    NK = FrontEdges(CF,LeftVertex)                  ! The corner of front edge called NK that an edge is defined upon it 
				    NV = FrontEdges(CF,RightVertex)				   ! The other corner of the current front edge called NV
				    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
				    ANL = Angles(CF,1)                              ! The Angle between current front edge and next front that share NK as common point 
				
	                print *,'ME:',ME			  
	                print *,'ANR1:',ANR*180/pi
	                print *,'ANL1:',ANL*180/pi

				    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,1,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Left)
				    Call DefineEdge(Dim,NC,NP,NK,NV,CF,ME,ANL,ANR,Corn,Neib,FrontEdges,Fronts,X,Y,States,newQuad,MFEC,NFE,NFEC,NF_Vertex_Left,1)
			
	                print *,'F:',CF,'&',':','newQuad1: ',newQuad(1),newQuad(2),newQuad(3),newQuad(4)
                
				    !------------------------------------- End of Definition of Side Edge on Right Vertex of Front Edge-------------------------------
	
				    Call TopEdgeRecovery(Dim,NC,newQuad(1),newQuad(2),Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
				
				    if(TopEdgeRecoveryIsPossible) then

					    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
                    
                        if(QFPossible) then
                        
					        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)   
					        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
					        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                        else
                    
                            print *,'---------->>>>>>>> Quad Formation Impossible! <<<<<<<<----------'
                                
                            Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                        Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                        
                        endif

                    else
					
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

				    endif
				
			    endif
					
		    endif

    ! ___________________________________________ End of Processing Fronts With Second Priority ___________________________________________

    ! ___________________________________________ Processing Fronts With Third Priority ____________________________________________

	    Case(STATUS_TWO)

	    !-------------------- Situation in which one side edge on the right should be defined ------------------------
!Part 8:
		    !------------------------------ Initializing Base Edge of new Quadrilateral ----------------------------------

		    newQuad(3)=FrontEdges(CF,LeftVertex)                   ! Saving point number of third corner of new Quad
		    newQuad(4)=FrontEdges(CF,RightVertex)                  ! Saving point number of fourth corner of new Quad

		    !------------ In this loop we find the corner of the Main Front that is not part of its Front Edge ------------

		    do J=1,3									   ! In this loop we find the corner of the Main Front that is not part of its Front Edge
			    if(Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,LeftVertex) .And. Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,RightVertex)) then
				    MFEC=Corn(FrontEdges(CF,Element),J)          ! MFEC(Main Front Element Corner)
				    MQE = Neib(FrontEdges(CF,Element),J)		 ! MQE(Main Quad Element of Current Front)
			    endif								   
		    end do

		    NK = FrontEdges(CF,LeftVertex)
		    NV = FrontEdges(CF,RightVertex)
		    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
		    ANL = Angles(CF,1)

		    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,1,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Left)

		    do J=1,3
			    if(Corn(NFE,J) == NFEC) then
				    NQE = Neib(NFE,J)
				    exit 
			    endif
		    end do

		    MFN = GetNorm(Dim,NK,NV,X,Y)                   ! Main Front Norm
		    NFN = GetNorm(Dim,NK,NF_Vertex_Left,X,Y) 	   ! Neibour Front Norm
		    Ratio = DMAX1(MFN,NFN)/DMIN1(MFN,NFN)

		    !----------------------------------- Finding Quad Element of Bigger Front ----------------------------------------
		    if(MFN > NFN) then
			    BQE = MQE
		    else
			    BQE = NQE
		    endif

		    Call GetSurroundingElements(Dim,Corn,Neib,NK,ME,TElms,QElms,TEC,QEC)

		    if((QEC>=5 .And. ANL<Alpha1) .Or. (QEC<5 .And. ANL<Alpha2)) then !--------------------------------------- Seaming -------------------------------------------

			    if(Ratio <= 2.5) then
   
				    print *,'------->>>>>>>>>>>> Seam Operation in STATUS_TWO <<<<<<<<<<<<-------'
!Part 9:				
				    Call Seam(Dim,NC,NP,NK,NV,NF_Vertex_Left,CF,NFI,Corn,Neib,X,Y,FrontEdges,Fronts,States,current_Level,2,SpecialCasePossible)
				
				    if(SpecialCasePossible) then
   
					    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
   
				    else
   
					    print *,'------------->>>>>>>>>>>>> Seam Impossible!!! <<<<<<<<<<<<<<<------------'
                    
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
   
				    endif

			    elseif(Ratio > 2.5) then

				    if(BQE/=0 .And. isNotChevron(Dim,Corn,X,Y,BQE)) then

					    print *,'------->>>>>>>>>>>> Transition Seam in STATUS_TWO <<<<<<<<<<<<-------'
!Part 10:                    
					    Call TransitionSeam(Dim,NC,NP,NBC,NK,NV,NF_Vertex_Left,CF,NFI,ME,NFE,MFEC,NFEC,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,SpecialCasePossible)
					
                        if(SpecialCasePossible) then
                        
                            Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                        else
                        
                            print *,'------->>>>>>>>>>>> Transition Seam is Impossible <<<<<<<<<<<<-------'
                        
                            Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                        Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                        
                        endif

                    else

                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

				    endif

			    endif

		    else

			    if(Ratio >= 2.5 .And. BQE/=0 .And. isNotChevron(Dim,Corn,X,Y,BQE)) then
				
				    print *,'------->>>>>>>>>>>> Transition Split in STATUS_TWO <<<<<<<<<<<<-------'
!Part 11:                
				    Call TransitionSplit(Dim,NC,NP,NBC,NK,NV,NF_Vertex_Left,ME,NFE,MFEC,NFEC,CF,NFI,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,SpecialCasePossible)
                
                    if(SpecialCasePossible) then
                    
				        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                
                    else
                
                        print *,'------->>>>>>>>>>>> Transition Split is Impossible <<<<<<<<<<<<-------'
                    
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                    endif

			    else
!Part 12:
				    newQuad(2)=NF_Vertex_Left								 ! Saving point number of second corner of new Quad

				    !------------------------------------- Defining Side Edge on Right Vertex of Front Edge ------------------------------------------
				
				    NK = FrontEdges(CF,RightVertex)                 ! The corner of front edge called NK that an edge is defined upon it 
				    NV = FrontEdges(CF,LeftVertex)				   ! The other corner of the current front edge called NV
				    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
				    ANR = Angles(CF,2)                              ! The Angle between current front edge and next front that share NK as common point
		
                    print *,'ME:',ME
		            print *,'ANR2:',ANR*180/pi
		            print *,'ANL2:',Angles(CF,1)*180/pi
        
				    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,2,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Right)

				    Call DefineEdge(Dim,NC,NP,NK,NV,CF,ME,ANL,ANR,Corn,Neib,FrontEdges,Fronts,X,Y,States,newQuad,MFEC,NFE,NFEC,NF_Vertex_Right,2)
				
		            print *,'F:',CF,'&',':','newQuad2: ',newQuad(1),newQuad(2),newQuad(3),newQuad(4)
                
				    !------------------------------------- End of Definition of Side Edge on Right Vertex of Front Edge-------------------------------

				    Call TopEdgeRecovery(Dim,NC,newQuad(1),newQuad(2),Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
				
				    if(TopEdgeRecoveryIsPossible) then

					    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
                    
                        if(QFPossible) then
                        
					        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
					        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
					        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                    
                        else
                    
                            print *,'---------->>>>>>>> Quad Formation Impossible! <<<<<<<<----------'
                                
                            Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                        Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                        
                        endif
	
                    else
					
                        Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                    Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                    Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

				    endif

			    endif

		    endif
    ! ___________________________________________ End of Processing Fronts With Third Priority ____________________________________________

    ! ____________________________________________ Processing Fronts With the Least Priority _________________________________________
	
	    Case(STATUS_ZERO)

	    !-------------------------- Situation in which two side edges should be defined ------------------------------
!Part 13:
		    !------------------------------ Initializing Base Edge of new Quadrilateral ----------------------------------

		    newQuad(3)=FrontEdges(CF,LeftVertex)                   ! Saving point number of third corner of new Quad
		    newQuad(4)=FrontEdges(CF,RightVertex)                  ! Saving point number of fourth corner of new Quad

		    !------------ In this loop we find the corner of the Main Front that is not part of its Front Edge ------------

		    do J=1,3				
			    if(Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,LeftVertex) .And. Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,RightVertex)) then
				    MFEC=Corn(FrontEdges(CF,Element),J)          ! MFEC(Main Front Element Corner)
			    endif								   
		    end do
			
		    !------------------------------------- Defining Side Edge on Right Vertex of Front Edge ------------------------------------------
	
		    NK = FrontEdges(CF,RightVertex)                 ! The corner of front edge called NK that an edge is defined upon it 
		    NV = FrontEdges(CF,LeftVertex)				   ! The other corner of the current front edge called NV
		    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
		    ANR = Angles(CF,2)                              ! The Angle between current front edge and next front that share NK as common point
		    print *,'ME:',ME
		    print *,'ANR0:',ANR*180/pi
		    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,2,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Right)
		    Call DefineEdge(Dim,NC,NP,NK,NV,CF,ME,ANL,ANR,Corn,Neib,FrontEdges,Fronts,X,Y,States,newQuad,MFEC,NFE,NFEC,NF_Vertex_Right,2)
		
		    !------------------------------------- End of Definition of Side Edge on Right Vertex of Front Edge-------------------------------
		
		    !------------ In this loop we find the corner of the Main Front that is not part of its Front Edge ------------

		    do J=1,3					
			    if(Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,LeftVertex) .And. Corn(FrontEdges(CF,Element),J) /= FrontEdges(CF,RightVertex)) then
				    MFEC=Corn(FrontEdges(CF,Element),J)          ! MFEC(Main Front Element Corner)
			    endif								   
		    end do

		    !------------------------------------- Defining Side Edge on Left Vertex of Front Edge--------------------------------------------

		    NK = FrontEdges(CF,LeftVertex)                  ! The corner of front edge called NK that an edge is defined upon it 
		    NV = FrontEdges(CF,RightVertex)				   ! The other corner of the current front edge called NV
		    ME = FrontEdges(CF,Element)					   ! Main Element of front edge
		    ANL = Angles(CF,1)                              ! The Angle between current front edge and next front that share NK as common point 
		    print *,'ANL0:',ANL*180/pi

		    Call getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,CF,1,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex_Left)
		    Call DefineEdge(Dim,NC,NP,NK,NV,CF,ME,ANL,ANR,Corn,Neib,FrontEdges,Fronts,X,Y,States,newQuad,MFEC,NFE,NFEC,NF_Vertex_Left,1)

		    !------------------------------------- End of Definition of Side Edge on Right Vertex of Front Edge-------------------------------
	
            print *,'F:',CF,'&',':','newQuad0: ',newQuad(1),newQuad(2),newQuad(3),newQuad(4)	
					
		    Call TopEdgeRecovery(Dim,NC,newQuad(1),newQuad(2),Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)

		    if(TopEdgeRecoveryIsPossible) then

			    Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
            
                if(QFPossible) then
                
			        Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
			        Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
			        Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
            
                else
            
                    print *,'---------->>>>>>>> Quad Formation Impossible! <<<<<<<<----------'
                                
                    Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	                Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	                Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	                Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	                Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)
                
                endif
		
            else
			
                Call MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,CF,Fronts,FrontEdges(CF,RightVertex),FrontEdges(CF,LeftVertex),newQuad)
	            Call QuadrilateralFormation(Dim,Corn,Neib,FrontEdges(CF,Element),newQuad,NC,QFPossible)
	            Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
	            Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
	            Call UpdateStates(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,Angles)

		    endif

	    !----------------------------------- End of situation in which two side edges should be defined --------------------------------

    !________________________________________________ End of Processing Fronts With the Least Priority ___________________________________

    End Select
    
endif
!===========================================================================================
End Subroutine DefineSideEdges
!*********************************************************************************************
