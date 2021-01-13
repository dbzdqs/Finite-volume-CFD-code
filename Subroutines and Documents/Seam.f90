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
Subroutine Seam(Dim,NC,NP,NK,NV,NF_Vertex,CF,NFI,Corn,Neib,X,Y,FrontEdges,Fronts,States,current_Level,Situation,SeamPossible)
Implicit None
!===========================================================================================
Intent(In)::Dim,NK,NV,NF_Vertex,CF,NFI,Situation,current_Level
Intent(InOut)::Corn,Neib,X,Y,FrontEdges,States,NC,NP,Fronts,SeamPossible

Integer,Parameter::to_Left = 1
Integer,Parameter::to_Right = 2
Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4
Integer,Parameter::PF=2	!---- PF(Processed Fronts)

Integer::Dim,Fronts,NC,NP,NK,NV,NF_Vertex,Nt,Nn,MFEC,CF,DTV,DTF,NFI,ME,MQ,FE,TEC_NV,QEC_NV,TEC_NF,QEC_NF,C,N1,N2,N3,N4,E,P1,P2,P3,index,I,J,It,Iv,Ik,current_Level,Situation,N,C_i,NV_i,NF_i,newN1,newN4,getNeibour
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::TElms_NV,QElms_NV,TEF_NV,QEF_NV,TElms_NF,QElms_NF,TEF_NF,QEF_NF
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::hasCorner,ElementInvertAt_NV,ElementInvertAt_NF,ElementInverted,SeamPossible,NumOfEdgesInTheLoopIsEven,isOnTheBoundary,LoopIsEvenToLeft,LoopIsEvenToRight,QFPossible
Real(8)::x_value,y_value,x_old_NV,y_old_NV,x_old_NF,y_old_NF
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
ElementInvertAt_NV = .False.
ElementInvertAt_NF = .False.

x_old_NV = X(NV)
y_old_NV = Y(NV)

x_old_NF = X(NF_Vertex)
y_old_NF = Y(NF_Vertex)

!Part 1:

!-------------------------- Recovering Edge Eo (Owen 1999) ---------------------------------

Call TopEdgeRecovery(Dim,NC,NV,NF_Vertex,Corn,Neib,X,Y,FrontEdges,Fronts,States,SeamPossible)

if(SeamPossible) then

!Part 2:
    
	!------------------------------- Determining Corner Nt -------------------------------------

	ME = FrontEdges(CF,Element)

	do I=1,3
		
		if(Corn(ME,I)/=NK .And. Corn(ME,I)/=NV) then
			
			MFEC = Corn(ME,I)
			MQ = Neib(ME,I)
			exit

		endif

	end do

	if(Corn(ME,1) == MFEC) then
		
		P1 = NF_Vertex 

	else

		P1 = Corn(ME,1)

	endif

	if(Corn(ME,2) == MFEC) then
		
		P2 = NF_Vertex 

	else

		P2 = Corn(ME,2)

	endif

	if(Corn(ME,3) == MFEC) then
		
		P3 = NF_Vertex 

	else

		P3 = Corn(ME,3)

	endif

	Call TriangleFormation(Dim,NC,Corn,Neib,ME,P1,P2,P3)

	ME = getNeibour(Dim,Corn,Neib,NK,NV,MQ)
	Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,ME)

	do I=1,3
		if(Corn(ME,I) == NK) then
			FE = Neib(ME,I) !------- FE(Facing Element to ME at NK)
			exit
		endif
	end do

	do I=1,3
		if(Neib(FE,I) == ME) then
			Nt = Corn(FE,I)
			exit
		endif
    end do

!Part 3:    
    
    Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,FE,Nn)

	do I=1,3
		if(Corn(ME,I) == NK) then
			FE = Neib(ME,I) !------- FE(Facing Element to ME at NK)
			exit
		endif
	end do

	do I=1,3
		if(Corn(FE,I) == NV) then
			DTF = Neib(FE,I) !--- DTF(Deleting Tri having corner NF_Vertex)
		endif
		if(Corn(FE,I) == NF_Vertex) then
			DTV = Neib(FE,I) !--- DTV(Deleting Tri having corner NV)
		endif
	end do
!Part 4:
	Select Case(Situation)
		
		Case(1)
			if(isOnTheBoundary(Dim,Fronts,Nt,FrontEdges,States)) then
				LoopIsEvenToLeft = NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,CF,NV,Nt,to_Left)
				LoopIsEvenToRight = NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,NFI,NF_Vertex,Nt,to_Right)
			endif
		Case(2)
			if(isOnTheBoundary(Dim,Fronts,Nt,FrontEdges,States)) then
				LoopIsEvenToRight = NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,CF,NV,Nt,to_Right)
				LoopIsEvenToLeft = NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,NFI,NF_Vertex,Nt,to_Left)
			endif
	End Select
!Part 5:
	!----------------------------- Composing the Quadrilateral ---------------------------------
	Select Case(Situation)

		Case(1)

			newQuad(1) = Nt
			newQuad(2) = NV
			newQuad(3) = NK
			newQuad(4) = NF_Vertex

		Case(2)

			newQuad(1) = Nt
			newQuad(2) = NF_Vertex
			newQuad(3) = NK
			newQuad(4) = NV 

	End Select

	Call QuadrilateralFormation(Dim,Corn,Neib,ME,newQuad,NC,QFPossible)

!Part 6:    
    
	!------------ Pre Calculating Surrounding Elements Faces of NV and NF_Vertex ---------------

	Call GetSurroundingElements(Dim,Corn,Neib,NV,NC,TElms_NV,QElms_NV,TEC_NV,QEC_NV)
	Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms_NV,QEC_NV,TElms_NV,TEC_NV,TEF_NV,QEF_NV)

	Call GetSurroundingElements(Dim,Corn,Neib,NF_Vertex,NC,TElms_NF,QElms_NF,TEC_NF,QEC_NF)
	Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms_NF,QEC_NF,TElms_NF,TEC_NF,TEF_NF,QEF_NF)

!Part 7:     
    
	!------------------------- Specifying Neibours of the Quadrilateral ----------------------------

	Select Case(Situation)
		
		Case(1)

			do I=1,4
				if(hasCorner(Dim,Corn,Neib(NC,I),Nt) .And. hasCorner(Dim,Corn,Neib(NC,I),NV)) then
					N1 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),NK) .And. hasCorner(Dim,Corn,Neib(NC,I),NV)) then
					N2 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),NK) .And. hasCorner(Dim,Corn,Neib(NC,I),NF_Vertex)) then
					N3 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),Nt) .And. hasCorner(Dim,Corn,Neib(NC,I),NF_Vertex)) then
					N4 = Neib(NC,I)
				endif 
			end do

		Case(2)

			do I=1,4
				if(hasCorner(Dim,Corn,Neib(NC,I),Nt) .And. hasCorner(Dim,Corn,Neib(NC,I),NF_Vertex)) then
					N1 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),NK) .And. hasCorner(Dim,Corn,Neib(NC,I),NF_Vertex)) then
					N2 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),NK) .And. hasCorner(Dim,Corn,Neib(NC,I),NV)) then
					N3 = Neib(NC,I)
				endif
				if(hasCorner(Dim,Corn,Neib(NC,I),Nt) .And. hasCorner(Dim,Corn,Neib(NC,I),NV)) then
					N4 = Neib(NC,I)
				endif 
			end do

        End Select

!Part 8:        
        
	!----------------------------------- Updating Fronts ---------------------------------------

	Call ManageFronts(Dim,Corn,Neib,NC,FrontEdges,Fronts,States,newQuad,current_Level)
    
!Part 9:    

	!------------ Special Case: N1 = N4 is a Quadrilateral adjacent at two sides -----------------
	do
		if(N1 /= N4) exit

		N = N1

		do I=1,4
			if(Corn(N,I)/=Nt .And. Corn(N,I)/=NV .And. Corn(N,I)/=NF_Vertex) then
				C = Corn(N,I)
				C_i = I
			endif
			if(Corn(N,I) == NV) NV_i = I
			if(Corn(N,I) == NF_Vertex) NF_i = I
		end do

		Select Case(Situation)

			Case(1)

				!------- Specifying new N1 ----------
				if(C_i > NV_i) then
					if(C_i - NV_i == 1) then
						newN1 = Neib(N,NV_i)
					else
						newN1 = Neib(N,C_i)
					endif
				else
					if(NV_i - C_i == 1) then
						newN1 = Neib(N,C_i)
					else
						newN1 = Neib(N,NV_i)
					endif
				endif
				!------- Specifying new N4 ----------
				if(C_i > NF_i) then
					if(C_i - NF_i == 1) then
						newN4 = Neib(N,NF_i)
					else
						newN4 = Neib(N,C_i)
					endif
				else
					if(NF_i - C_i == 1) then
						newN4 = Neib(N,C_i)
					else
						newN4 = Neib(N,NF_i)
					endif
				endif

			Case(2)

				!------- Specifying new N1 ----------
				if(C_i > NF_i) then
					if(C_i - NF_i == 1) then
						newN1 = Neib(N,NF_i)
					else
						newN1 = Neib(N,C_i)
					endif
				else
					if(NF_i - C_i == 1) then
						newN1 = Neib(N,C_i)
					else
						newN1 = Neib(N,NF_i)
					endif
				endif
				!------- Specifying new N4 ----------
				if(C_i > NV_i) then
					if(C_i - NV_i == 1) then
						newN4 = Neib(N,NV_i)
					else
						newN4 = Neib(N,C_i)
					endif
				else
					if(NV_i - C_i == 1) then
						newN4 = Neib(N,C_i)
					else
						newN4 = Neib(N,NV_i)
					endif
				endif

		End Select

		!------------------ Modifying NC neibours ------------------
		do I=1,4
			if(Neib(NC,I) == N1) Neib(NC,I) = newN1
			if(Neib(NC,I) == N4) Neib(NC,I) = newN4
		end do

		N1 = newN1
		N4 = newN4

		!-------------- Modify new N1 and new N4 Neibours ----------
		do I=1,4
			if(Neib(newN1,I) == N) Neib(newN1,I) = NC
			if(Neib(newN4,I) == N) Neib(newN4,I) = NC
		end do

		!------------------------ Deleting N -----------------------
		do I=1,4
			Corn(N,I) = Corn(N,1)
		end do

		Nt = C

    end do

!Part 10:    
    
	!---------------------------- Merging two Corners at midway --------------------------------

	x_value = (X(NF_Vertex) + X(NV))/2
	y_value = (Y(NF_Vertex) + Y(NV))/2

	X(NF_Vertex) = x_value
	Y(NF_Vertex) = y_value
	X(NV) = x_value
	Y(NV) = y_value

!Part 11:    
    
	!----------------- Checking Inversion of Elements Surrounding NV and NF_Vertex -----------------

	ElementInvertAt_NV = ElementInverted(Dim,Corn,X,Y,TElms_NV,QElms_NV,TEC_NV,QEC_NV,TEF_NV,QEF_NV)

	if(ElementInvertAt_NV) then

		X(NV) = x_old_NV
		Y(NV) = y_old_NV
		X(NF_Vertex) = x_old_NV
		Y(NF_Vertex) = y_old_NV

	endif

	ElementInvertAt_NF = ElementInverted(Dim,Corn,X,Y,TElms_NF,QElms_NF,TEC_NF,QEC_NF,TEF_NF,QEF_NF) 

	if(ElementInvertAt_NF) then
		
		X(NF_Vertex) = x_old_NF
		Y(NF_Vertex) = y_old_NF
		X(NV) = x_old_NF
		Y(NV) = y_old_NF

	endif

	!----------------------------- Final Inversion Checking -------------------------------------

	ElementInvertAt_NV = ElementInverted(Dim,Corn,X,Y,TElms_NV,QElms_NV,TEC_NV,QEC_NV,TEF_NV,QEF_NV)
	ElementInvertAt_NF = ElementInverted(Dim,Corn,X,Y,TElms_NF,QElms_NF,TEC_NF,QEC_NF,TEF_NF,QEF_NF)

	if(ElementInvertAt_NV .Or. ElementInvertAt_NF) then
		
		X(NV) = x_old_NV
		Y(NV) = y_old_NV
		X(NF_Vertex) = x_old_NF
		Y(NF_Vertex) = y_old_NF

		if((.Not. LoopIsEvenToLeft) .Or. (.Not. LoopIsEvenToRight)) then

			!------------------- Modifying Nt corner to Nn -------------	
			do I=1,4
				if(Corn(NC,I) == Nt) then
					Corn(NC,I) = Nn
					exit
				endif
			end do

			!---------------------- Modifying Neibours --------------
			Select Case(Situation)
				
				Case(1)
					
					!---------- Modify NC Neibours ----------
					do I=1,4
						if(Neib(NC,I) == N1) Neib(NC,I) = DTV
						if(Neib(NC,I) == N4) Neib(NC,I) = DTF
					end do

					!--------- Modify N1,N4 Neibours --------
					do I=1,4
						if(Neib(N1,I) == NC) Neib(N1,I) = DTV
						if(Neib(N4,I) == NC) Neib(N4,I) = DTF
					end do
					
					!------- Reconstructing DTV,DTF ---------
					do I=1,3
						if(Neib(DTV,I) == FE) Corn(DTV,I) = Nt
						if(Neib(DTV,I) == DTF) Corn(DTV,I) = NV
						if(Neib(DTV,I) == N1) Corn(DTV,I) = Nn

						if(Neib(DTF,I) == FE) Corn(DTF,I) = Nt
						if(Neib(DTF,I) == DTV) Corn(DTF,I) = NF_Vertex
						if(Neib(DTF,I) == N4) Corn(DTF,I) = Nn
					end do
					!------- Modify DTV,DTF Neibours --------
					do I=1,3
						if(Neib(DTV,I) == FE) Neib(DTV,I) = NC
						if(Neib(DTF,I) == FE) Neib(DTF,I) = NC
					end do

					!------------------ Update newly added Fronts ----------------
					do I=1,Fronts
						if(States(I) /= Processed) then
							if(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NF_Vertex) then
								FrontEdges(I,LeftVertex) = Nn
								FrontEdges(I,Element) = DTF
								print*,'UP1: ',I,FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex) 
							endif

							if(FrontEdges(I,LeftVertex)==NV .And. FrontEdges(I,RightVertex)==Nt) then
								FrontEdges(I,RightVertex) = Nn
								FrontEdges(I,Element) = DTV
								print*,'UP2: ',I,FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex) 
							endif
						endif
					end do

				Case(2)

					!---------- Modify NC Neibours ----------
					do I=1,4
						if(Neib(NC,I) == N1) Neib(NC,I) = DTF
						if(Neib(NC,I) == N4) Neib(NC,I) = DTV
					end do

					!--------- Modify N1,N4 Neibours --------
					do I=1,4
						if(Neib(N1,I) == NC) Neib(N1,I) = DTF
						if(Neib(N4,I) == NC) Neib(N4,I) = DTV
					end do
					
					!------- Reconstructing DTV,DTF ---------
					do I=1,3
						if(Neib(DTV,I) == FE) Corn(DTV,I) = Nt
						if(Neib(DTV,I) == DTF) Corn(DTV,I) = NV
						if(Neib(DTV,I) == N4) Corn(DTV,I) = Nn

						if(Neib(DTF,I) == FE) Corn(DTF,I) = Nt
						if(Neib(DTF,I) == DTV) Corn(DTF,I) = NF_Vertex
						if(Neib(DTF,I) == N1) Corn(DTF,I) = Nn
					end do
					!------- Modify DTV,DTF Neibours --------
					do I=1,3
						if(Neib(DTV,I) == FE) Neib(DTV,I) = NC
						if(Neib(DTF,I) == FE) Neib(DTF,I) = NC
					end do

					!------------------ Update newly added Fronts ----------------
					do I=1,Fronts
						if(States(I) /= Processed) then
							if(FrontEdges(I,LeftVertex)==NF_Vertex .And. FrontEdges(I,RightVertex)==Nt) then
								FrontEdges(I,RightVertex) = Nn
								FrontEdges(I,Element) = DTF
								print*,'UP1: ',I,FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex) 
							endif

							if(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NV) then
								FrontEdges(I,LeftVertex) = Nn
								FrontEdges(I,Element) = DTV
								print*,'UP2: ',I,FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex) 
							endif
						endif
					end do

			End Select

		endif

	else

		Select Case(Situation)

			Case (1)

				do I=1,Fronts

					if(FrontEdges(I,RightVertex) == NV) then
						index = I	
					endif

					if(States(I) /= Processed) then

						if(Corn(N1,4) == 0) then
							if(FrontEdges(I,LeftVertex) == NV .And. FrontEdges(I,RightVertex) == Nt) then
								print*,'Seam Front Processed: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								States(I) = Processed
							endif
						endif

						if(Corn(N4,4) == 0) then
							if(FrontEdges(I,LeftVertex) == Nt .And. FrontEdges(I,RightVertex) == NF_Vertex) then
								print*,'Seam Front Processed: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								States(I) = Processed
							endif
						endif

					endif

				end do

				if(Corn(N1,4) /= 0 .And. Corn(N4,4) == 0) then
					do I=1,Fronts
						if(States(I) == Processed) then
							if(FrontEdges(I,LeftVertex)==NV .And. FrontEdges(I,RightVertex)==Nt) then
								States(I) = -2
								FrontEdges(I,Element) = N4
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							elseif(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NV) then
								States(I) = -2
								FrontEdges(I,Element) = N4
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							endif
						endif
					end do
				endif

				if(Corn(N4,4) /= 0 .And. Corn(N1,4) == 0) then
					do I=1,Fronts
						if(States(I) == Processed) then
							if(FrontEdges(I,LeftVertex)==NF_Vertex .And. FrontEdges(I,RightVertex)==Nt) then
								States(I) = -2
								FrontEdges(I,Element) = N1
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							elseif(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NF_Vertex) then
								States(I) = -2
								FrontEdges(I,Element) = N1
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							endif
						endif
					end do
				endif

				FrontEdges(index,RightVertex) = NF_Vertex

			Case (2)

				do I=1,Fronts

					if(FrontEdges(I,LeftVertex) == NV) then
						index = I	
					endif

					if(States(I) /= Processed) then

						if(Corn(N1,4) == 0) then
							if(FrontEdges(I,LeftVertex) == NF_Vertex .And. FrontEdges(I,RightVertex) == Nt) then
								print*,'Seam Front Processed: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								States(I) = Processed
							endif
						endif

						if(Corn(N4,4) == 0) then
							if(FrontEdges(I,LeftVertex) == Nt .And. FrontEdges(I,RightVertex) == NV) then
								print*,'Seam Front Processed: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								States(I) = Processed
							endif
						endif

					endif

				end do

				if(Corn(N1,4) /= 0 .And. Corn(N4,4) == 0) then
					do I=1,Fronts
						if(States(I) == Processed) then
							if(FrontEdges(I,LeftVertex)==NF_Vertex .And. FrontEdges(I,RightVertex)==Nt) then
								States(I) = -2
								FrontEdges(I,Element) = N4
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							elseif(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NF_Vertex) then
								States(I) = -2
								FrontEdges(I,Element) = N4
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							endif
						endif
					end do
				endif

				if(Corn(N4,4) /= 0 .And. Corn(N1,4) == 0) then
					do I=1,Fronts
						if(States(I) == Processed) then
							if(FrontEdges(I,LeftVertex)==NV .And. FrontEdges(I,RightVertex)==Nt) then
								States(I) = -2
								FrontEdges(I,Element) = N1
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							elseif(FrontEdges(I,LeftVertex)==Nt .And. FrontEdges(I,RightVertex)==NV) then
								States(I) = -2
								FrontEdges(I,Element) = N1
								print*,'Seam Front Added: ',I,' : ',FrontEdges(I,LeftVertex),FrontEdges(I,RightVertex)
								exit
							endif
						endif
					end do
				endif

				FrontEdges(index,LeftVertex) = NF_Vertex

		End Select

		Select Case(Situation)
			
			Case(1)
				
				!-------------- Updating N1 Neibours -------
				if(Corn(N1,4) /= 0) then
					do I=1,4
						if(Corn(N1,I) == NV) Iv = I
						if(Corn(N1,I) == Nt) It = I
					end do

					if(Iv > It) then
						if(Iv - It == 1) then
							Neib(N1,It) = N4
						else
							Neib(N1,Iv) = N4
						endif
					else
						if(It - Iv == 1) then
							Neib(N1,Iv) = N4
						else
							Neib(N1,It) = N4
						endif
					endif
				else
					do I=1,3
						if(Corn(N1,I)/=NV .And. Corn(N1,I)/=Nt) then
							Neib(N1,I) = N4
							exit
						endif
					end do
				endif
				!-------------- Updating N2 Neibours -------
				if(Corn(N2,4) /= 0) then
					do I=1,4
						if(Corn(N2,I) == NV) Iv = I
						if(Corn(N2,I) == NK) Ik = I
					end do

					if(Iv > Ik) then
						if(Iv - Ik == 1) then
							Neib(N2,Ik) = N3
						else
							Neib(N2,Iv) = N3
						endif
					else
						if(Ik - Iv == 1) then
							Neib(N2,Iv) = N3
						else
							Neib(N2,Ik) = N3
						endif
					endif
				else
					do I=1,3
						if(Corn(N2,I)/=NV .And. Corn(N2,I)/=NK) then
							Neib(N2,I) = N3
							exit
						endif
					end do
				endif
				!-------------- Updating N3 Neibours -------
				if(Corn(N3,4) /= 0) then
					do I=1,4
						if(Corn(N3,I) == NF_Vertex) Iv = I
						if(Corn(N3,I) == NK) Ik = I
					end do

					if(Iv > Ik) then
						if(Iv - Ik == 1) then
							Neib(N3,Ik) = N2
						else
							Neib(N3,Iv) = N2
						endif
					else
						if(Ik - Iv == 1) then
							Neib(N3,Iv) = N2
						else
							Neib(N3,Ik) = N2
						endif
					endif
				else
					do I=1,3
						if(Corn(N3,I)/=NF_Vertex .And. Corn(N3,I)/=NK) then
							Neib(N3,I) = N2
							exit
						endif
					end do
				endif
				!-------------- Updating N4 Neibours -------
				if(Corn(N4,4) /= 0) then
					do I=1,4
						if(Corn(N4,I) == NF_Vertex) Iv = I
						if(Corn(N4,I) == Nt) It = I
					end do

					if(Iv > It) then
						if(Iv - It == 1) then
							Neib(N4,It) = N1
						else
							Neib(N4,Iv) = N1
						endif
					else
						if(It - Iv == 1) then
							Neib(N4,Iv) = N1
						else
							Neib(N4,It) = N1
						endif
					endif
				else
					do I=1,3
						if(Corn(N4,I)/=NF_Vertex .And. Corn(N4,I)/=Nt) then
							Neib(N4,I) = N1
							exit
						endif
					end do
				endif

			Case(2)

				!-------------- Updating N1 Neibours -------
				if(Corn(N1,4) /= 0) then
					do I=1,4
						if(Corn(N1,I) == NF_Vertex) Iv = I
						if(Corn(N1,I) == Nt) It = I
					end do

					if(Iv > It) then
						if(Iv - It == 1) then
							Neib(N1,It) = N4
						else
							Neib(N1,Iv) = N4
						endif
					else
						if(It - Iv == 1) then
							Neib(N1,Iv) = N4
						else
							Neib(N1,It) = N4
						endif
					endif
				else
					do I=1,3
						if(Corn(N1,I)/=NF_Vertex .And. Corn(N1,I)/=Nt) then
							Neib(N1,I) = N4
							exit
						endif
					end do
				endif
				!-------------- Updating N2 Neibours -------
				if(Corn(N2,4) /= 0) then
					do I=1,4
						if(Corn(N2,I) == NF_Vertex) Iv = I
						if(Corn(N2,I) == NK) Ik = I
					end do

					if(Iv > Ik) then
						if(Iv - Ik == 1) then
							Neib(N2,Ik) = N3
						else
							Neib(N2,Iv) = N3
						endif
					else
						if(Ik - Iv == 1) then
							Neib(N2,Iv) = N3
						else
							Neib(N2,Ik) = N3
						endif
					endif
				else
					do I=1,3
						if(Corn(N2,I)/=NF_Vertex .And. Corn(N2,I)/=NK) then
							Neib(N2,I) = N3
							exit
						endif
					end do
				endif
				!-------------- Updating N3 Neibours -------
				if(Corn(N3,4) /= 0) then
					do I=1,4
						if(Corn(N3,I) == NV) Iv = I
						if(Corn(N3,I) == NK) Ik = I
					end do

					if(Iv > Ik) then
						if(Iv - Ik == 1) then
							Neib(N3,Ik) = N2
						else
							Neib(N3,Iv) = N2
						endif
					else
						if(Ik - Iv == 1) then
							Neib(N3,Iv) = N2
						else
							Neib(N3,Ik) = N2
						endif
					endif
				else
					do I=1,3
						if(Corn(N3,I)/=NV .And. Corn(N3,I)/=NK) then
							Neib(N3,I) = N2
							exit
						endif
					end do
				endif
				!-------------- Updating N4 Neibours -------
				if(Corn(N4,4) /= 0) then
					do I=1,4
						if(Corn(N4,I) == NV) Iv = I
						if(Corn(N4,I) == Nt) It = I
					end do

					if(Iv > It) then
						if(Iv - It == 1) then
							Neib(N4,It) = N1
						else
							Neib(N4,Iv) = N1
						endif
					else
						if(It - Iv == 1) then
							Neib(N4,Iv) = N1
						else
							Neib(N4,It) = N1
						endif
					endif
				else
					do I=1,3
						if(Corn(N4,I)/=NV .And. Corn(N4,I)/=Nt) then
							Neib(N4,I) = N1
							exit
						endif
					end do
				endif

		End Select
		!*+*+*+*+*+

		do I=1,4
			Corn(NC,I) = NF_Vertex
		end do

		do I=1,TEC_NV
			E = TElms_NV(I)
			do J=1,3
				if(Corn(E,J) == NV) then
					Corn(E,J) = NF_Vertex
				endif  
			end do
		end do
		do I=1,QEC_NV
			E = QElms_NV(I)
			do J=1,4
				if(Corn(E,J) == NV) then
					Corn(E,J) = NF_Vertex
				endif  
			end do
		end do
!>>>>>>>>>>>>>>>>>>>>>>>>>>>+++++++++++++++++++++++++++++++++++++++++++++++
		do I=1,Fronts
			if(States(I) /= Processed) then
				if(FrontEdges(I,LeftVertex) == NV) then
					FrontEdges(I,LeftVertex) = NF_Vertex 
				endif

				if(FrontEdges(I,RightVertex) == NV) then
					FrontEdges(I,RightVertex) = NF_Vertex
				endif
			endif
		end do
!>>>>>>>>>>>>>>>>>>>>>>>>>>>+++++++++++++++++++++++++++++++++++++++++++++++
	endif

endif
!===========================================================================================
End Subroutine Seam 
!*********************************************************************************************
