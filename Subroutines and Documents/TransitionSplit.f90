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
Subroutine TransitionSplit(Dim,NC,NP,NBC,NK,NV,NF_Vertex,ME,NFE,MFEC,NFEC,CF,NFI,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,TopEdgeRecoveryIsPossible)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBC,NK,NV,NF_Vertex,ME,NFE,MFEC,NFEC,BFP,MFN,NFN,CF,NFI
Intent(InOut)::Corn,Neib,FrontEdges,X,Y,States,Fronts,TopEdgeRecoveryIsPossible

Integer,Parameter::Recently_Added = -2
Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4
Integer,Parameter::PF=2	!---- PF(Processed Fronts)
Integer,Parameter::FC=1

Integer::Dim,NC,NP,NBC,NK,NV,NF_Vertex,ME,NFE,MFEC,NFEC,CF,NFI,Fronts,N1,N2,I,I1,I2,J1,J2,J3,Q,newQ,newT1,newT2,C,C_i,NV_i,NF_Vertex_i,NK_i,P_i,E1,E2,getNeibour
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::TopEdgeRecoveryIsPossible,QFPossible
Real(8)::MFN,NFN,x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!____________________________ Case when ME and its adjacent Quad are split _________________

if(MFN > NFN) then
!Part 1:    
	print*,'>>>> ME and its adjacent Quad are split <<<<'
	!-------------------------- Finding the Quad (the Splitting Quad) ---------------------- 
	do I=1,3
		if(Corn(ME,I) == MFEC) then
			Q = Neib(ME,I)
			exit
		endif
	end do

	!----------------------- Finding the other Adjacent Corner to NV in Q ------------------
	do I=1,4
		if(Corn(Q,I) == NV) then
			NV_i = I
			if(I==1) then
				I1 = 4
			else
				I1 = I - 1
			endif

			if(I==4) then
				I2 = 1
			else
				I2 = I + 1
			endif
			exit
		endif
	end do

	if(Corn(Q,I1) /= NK) then
		C = Corn(Q,I1)
		C_i = I1
		NK_i = I2
	else
		C = Corn(Q,I2)
		C_i = I2
		NK_i = I1
	endif

	!------------------- Specifying Index of Other Corner of Splitting Quad ----------------
	
	do I=1,4
		if(I/=NV_i .And. I/=NK_i .And. I/=C_i) then
			P_i = I
			exit
		endif
	end do

	!-------------------------------- Specifying Indices of ME -----------------------------

	do I=1,3
		if(Corn(ME,I) == MFEC) J1 = I !---------- Coresponding Index of MFEC ---------------
		if(Corn(ME,I) == NV) J2 = I   !---------- Coresponding Index of NV -----------------
		if(Corn(ME,I) == NK) J3 = I   !---------- Coresponding Index of NK -----------------
	end do
!Part 2:
	!--------------------------- N1: Mid Point of Larger Front -----------------------------
	NP = NP + 1
	N1 = NP
	X(N1) = (X(NK)+X(NV))/2
	Y(N1) = (Y(NK)+Y(NV))/2

	!--------------------------- N2: Centroid Point of Bigger Quad -------------------------  
	NP = NP + 1
	N2 = NP
	Call CalcCentroidOfQuad(Dim,Corn,X,Y,Q,x_value,y_value)
	X(N2) = x_value	
	Y(N2) = y_value
	!------------------------------- Composing New Quad ------------------------------------
	NC = NC + 1
	newQ = NC
	Corn(newQ,NK_i) = N1
	Corn(newQ,NV_i) = NV
	Corn(newQ,C_i) = C
	Corn(newQ,P_i) = N2
	!--------------------- Composing new Triangle newT1(N2,NK,N1) --------------------------
	NC = NC + 1
	newT1 = NC
	Corn(newT1,J1) = N2
	Corn(newT1,J2) = NK
	Corn(newT1,J3) = N1
	Corn(newT1,4) = 0
	!--------------------- Composing new Triangle newT2(N1,NK,MFEC) ------------------------
	NC = NC + 1
	newT2 = NC
	Corn(newT2,J1) = MFEC
	Corn(newT2,J2) = N1
	Corn(newT2,J3) = NK
	Corn(newT2,4) = 0
	!----------------------------- Determining newQ neibours -------------------------------
    
    Call setNeibour(Dim,Corn,Neib,newQ,NV,N1,ME)
	E1 = getNeibour(Dim,Corn,Neib,NV,C,Q)
    Call setNeibour(Dim,Corn,Neib,newQ,NV,C,E1)
	Call setNeibour(Dim,Corn,Neib,newQ,N2,C,Q)
	Call setNeibour(Dim,Corn,Neib,newQ,N1,N2,newT1) 
	
	!---------------------- Updating Neibours of Neibours of newQ --------------------------

    Call setNeibour(Dim,Corn,Neib,E1,NV,C,newQ)

	!----------------------------- Determining newT1 neibours ------------------------------
	Call setNeibour(Dim,Corn,Neib,newT1,NK,N1,newT2)
    Call setNeibour(Dim,Corn,Neib,newT1,N1,N2,newQ)
	Call setNeibour(Dim,Corn,Neib,newT1,NK,N2,Q)
	Neib(newT1,4) = 0 
	!----------------------------- Determining newT2 neibours ------------------------------
	Call setNeibour(Dim,Corn,Neib,newT2,NK,N1,newT1)
    E2 = getNeibour(Dim,Corn,Neib,NK,MFEC,ME) 
    Call setNeibour(Dim,Corn,Neib,newT2,NK,MFEC,E2)
	Call setNeibour(Dim,Corn,Neib,newT2,N1,MFEC,ME)
	Neib(newT2,4) = 0
	!---------------------- Updating Neibours of Neibours of newT2 -------------------------
    
	Call setNeibour(Dim,Corn,Neib,E2,NK,MFEC,newT2)
!Part 3:    
	!--------------------------- Modifying the Splitting Quad ------------------------------

	Corn(Q,NV_i) = N2

    Call setNeibour(Dim,Corn,Neib,Q,NK,N2,newT1)
    Call setNeibour(Dim,Corn,Neib,Q,C,N2,newQ)
			    
	!------------------------------- Modifying ME triangle ---------------------------------
	do I=1,3
		if(Corn(ME,I) == NK) then 
            Corn(ME,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,ME,N1,MFEC,newT2)
    Call setNeibour(Dim,Corn,Neib,ME,N1,NV,newQ)
	!-------------------- Composing the newQ2(NF_vertex,N1,N2,NK) Quad --------------------- 
!Part 4:    
	Call TopEdgeRecovery(Dim,NC,NF_vertex,N1,Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
	if(TopEdgeRecoveryIsPossible) then
!Part 5:
		newQuad(NK_i) = NF_vertex
		newQuad(NV_i) = N1
		newQuad(C_i) = N2
		newQuad(P_i) = NK
		Call QuadrilateralFormation(Dim,Corn,Neib,newT1,newQuad,NC,QFPossible)		
		Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
!Part 6:
        !------------------------------- Modify split Front Info -------------------------------

	    print *,'Smaller Front Removed!!! ',NFI,': ',FrontEdges(NFI,LeftVertex),FrontEdges(NFI,RightVertex)
	    States(NFI) = Processed

	    if(FrontEdges(CF,LeftVertex) == NK) then

		    FrontEdges(CF,LeftVertex) = N1
		    Fronts = Fronts + 1
		    FrontEdges(Fronts,LeftVertex) = NF_vertex
		    FrontEdges(Fronts,RightVertex) = N1
		    FrontEdges(Fronts,Level) = FrontEdges(CF,Level)
		    States(Fronts) = Recently_Added
		    print *,'Front Added!!! ',Fronts,': ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)
		
	    else

		    FrontEdges(CF,RightVertex) = N1
		    Fronts = Fronts + 1
		    FrontEdges(Fronts,LeftVertex) = N1
		    FrontEdges(Fronts,RightVertex) = NF_vertex
		    FrontEdges(Fronts,Level) = FrontEdges(CF,Level)
		    States(Fronts) = Recently_Added
		    print *,'Front Added!!! ',Fronts,': ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

        endif
    
        !---------------------------------- Updating CF element --------------------------------
        
        FrontEdges(CF,Element) = getNeibour(Dim,Corn,Neib,NV,N1,newQ) 
        
		!--------------------------- Specifying new Front's Element ---------------------------
		do I=1,4
			if(Corn(NC,I) == NF_vertex) I1 = I
			if(Corn(NC,I) == N1) I2 = I
		end do

		if(I1 < I2) then
			if(I2 - I1 == 1) then
				FrontEdges(Fronts,Element) = Neib(NC,I1)
			else
				FrontEdges(Fronts,Element) = Neib(NC,I2)
			endif
		else
			if(I1 - I2 == 1) then
				FrontEdges(Fronts,Element) = Neib(NC,I2)
			else
				FrontEdges(Fronts,Element) = Neib(NC,I1)
			endif
        endif
        
        !--------------------------------- Smoothing Q -----------------------------------------
	    newQuad(NK_i) = Corn(Q,NK_i)
	    newQuad(NV_i) = Corn(Q,NV_i)
	    newQuad(C_i) = Corn(Q,C_i)
	    newQuad(P_i) = Corn(Q,P_i)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Q,X,Y)

	    !------------------------------- Smoothing newQ ----------------------------------------
	    newQuad(NK_i) = Corn(newQ,NK_i)
	    newQuad(NV_i) = Corn(newQ,NV_i)
	    newQuad(C_i) = Corn(newQ,C_i)
	    newQuad(P_i) = Corn(newQ,P_i)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,newQ,X,Y)
        
    else !-------------------- Case when edge recovery is not possible ----------------------
!Part 7:    
        !-------------------------- Reconstructing Q element --------------------------------
        do I=1,4
            if(Corn(Q,I) == N2) then
                Corn(Q,I) = NV
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,Q,C,NV,E1)
        Call setNeibour(Dim,Corn,Neib,Q,NK,NV,ME)
        !-------------------------- Reconstructing ME element -------------------------------
        do I=1,3
            if(Corn(ME,I) == N1) then
                Corn(ME,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,ME,NK,NV,Q)
        Call setNeibour(Dim,Corn,Neib,ME,NK,MFEC,E2)
        
        !-------------------------- Modifying E1,E2 neibours -------------------------------
        
        Call setNeibour(Dim,Corn,Neib,E1,C,NV,Q)
        Call setNeibour(Dim,Corn,Neib,E2,NK,MFEC,ME)
        
        !-------------------------- Deleting newT1,newT2,newQ ------------------------------
        
        Call DeleteCell(Dim,newT1,Corn)
        Call DeleteCell(Dim,newT2,Corn)
        Call DeleteCell(Dim,newQ,Corn)
        
	endif

!___________________________ Case when NFE and its adjacent Quad are split _________________	 		
else
!Part 1:    
	print*,'>>>> NFE and its adjacent Quad are split <<<<'
	!-------------------------- Finding the Quad (the Splitting Quad) ----------------------
	do I=1,3
		if(Corn(NFE,I) == NFEC) then
			Q = Neib(NFE,I)
			exit
		endif
	end do

	!------------------ Finding the other Adjacent Corner to NF_Vertex in Q ----------------
	do I=1,4
		if(Corn(Q,I) == NF_Vertex) then
			NF_Vertex_i = I
			if(I==1) then
				I1 = 4
			else
				I1 = I - 1
			endif

			if(I==4) then
				I2 = 1
			else
				I2 = I + 1
			endif
			exit
		endif
	end do

	if(Corn(Q,I1) /= NK) then
		C = Corn(Q,I1)
		C_i = I1
		NK_i = I2
	else
		C = Corn(Q,I2)
		C_i = I2
		NK_i = I1
	endif
	
	!------------------- Specifying Index of Other Corner of Splitting Quad ----------------
	
	do I=1,4
		if(I/=NF_Vertex_i .And. I/=NK_i .And. I/=C_i) then
			P_i = I
			exit
		endif
	end do

	!-------------------------------- Specifying Indices of ME -----------------------------

	do I=1,3
		if(Corn(NFE,I) == NFEC) J1 = I	    !-------- Coresponding Index of NFEC -----------
		if(Corn(NFE,I) == NF_Vertex) J2 = I !-------- Coresponding Index of NF_Vertex ------
		if(Corn(NFE,I) == NK) J3 = I        !-------- Coresponding Index of NK -------------
	end do
!Part 2:
	!--------------------------- N1: Mid Point of Larger Front -----------------------------
	NP = NP + 1
	N1 = NP
	X(N1) = (X(NK)+X(NF_Vertex))/2
	Y(N1) = (Y(NK)+Y(NF_Vertex))/2

	!--------------------------- N2: Centroid Point of Bigger Quad -------------------------  
	NP = NP + 1
	N2 = NP
	Call CalcCentroidOfQuad(Dim,Corn,X,Y,Q,x_value,y_value)
	X(N2) = x_value	
	Y(N2) = y_value
	!------------------------------- Composing New Quad ------------------------------------
	NC = NC + 1
	newQ = NC
	Corn(newQ,NK_i) = N1
	Corn(newQ,NF_Vertex_i) = NF_Vertex
	Corn(newQ,C_i) = C
	Corn(newQ,P_i) = N2
	!--------------------- Composing new Triangle newT1(N2,NK,N1) --------------------------
	NC = NC + 1
	newT1 = NC
	Corn(newT1,J1) = N2
	Corn(newT1,J2) = NK
	Corn(newT1,J3) = N1
	Corn(newT1,4) = 0
	!--------------------- Composing new Triangle newT2(N1,NK,MFEC) ------------------------
	NC = NC + 1
	newT2 = NC
	Corn(newT2,J1) = NFEC
	Corn(newT2,J2) = N1
	Corn(newT2,J3) = NK
	Corn(newT2,4) = 0
	!----------------------------- Determining newQ neibours -------------------------------
    
    Call setNeibour(Dim,Corn,Neib,newQ,NF_Vertex,N1,NFE)
    E1 = getNeibour(Dim,Corn,Neib,NF_Vertex,C,Q)
    Call setNeibour(Dim,Corn,Neib,newQ,NF_Vertex,C,E1)
	Call setNeibour(Dim,Corn,Neib,newQ,N2,C,Q)
    Call setNeibour(Dim,Corn,Neib,newQ,N1,N2,newT1)
	
	!---------------------- Updating Neibours of Neibours of newQ --------------------------

    Call setNeibour(Dim,Corn,Neib,E1,NF_Vertex,C,newQ)

	!----------------------------- Determining newT1 neibours ------------------------------
	Call setNeibour(Dim,Corn,Neib,newT1,N1,NK,newT2)
    Call setNeibour(Dim,Corn,Neib,newT1,N1,N2,newQ)
    Call setNeibour(Dim,Corn,Neib,newT1,NK,N2,Q)
	Neib(newT1,4) = 0 
	!----------------------------- Determining newT2 neibours ------------------------------
	Call setNeibour(Dim,Corn,Neib,newT2,N1,NK,newT1)
    E2 = getNeibour(Dim,Corn,Neib,NK,NFEC,NFE)
    Call setNeibour(Dim,Corn,Neib,newT2,NK,NFEC,E2)
	Call setNeibour(Dim,Corn,Neib,newT2,N1,NFEC,NFE)
	Neib(newT2,4) = 0
	!---------------------- Updating Neibours of Neibours of newT2 -------------------------
	
    Call setNeibour(Dim,Corn,Neib,E2,NK,NFEC,newT2)
!Part 3:
	!--------------------------- Modifying the Splitting Quad ------------------------------

	Corn(Q,NF_Vertex_i) = N2

    Call setNeibour(Dim,Corn,Neib,Q,NK,N2,newT1)
    Call setNeibour(Dim,Corn,Neib,Q,C,N2,newQ)
		    
	!------------------------------- Modifying NFE triangle ---------------------------------
	do I=1,3
		if(Corn(NFE,I) == NK) then 
            Corn(NFE,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,NFE,N1,NFEC,newT2)
    Call setNeibour(Dim,Corn,Neib,NFE,N1,NF_Vertex,newQ)

	!-------------------- Composing the newQ2(NF_vertex,N1,N2,NK) Quad ---------------------
!Part 4:
	Call TopEdgeRecovery(Dim,NC,NV,N1,Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
	if(TopEdgeRecoveryIsPossible) then
!Part 5:         
		newQuad(NK_i) = NV
		newQuad(NF_Vertex_i) = N1
		newQuad(C_i) = N2
		newQuad(P_i) = NK       
		Call QuadrilateralFormation(Dim,Corn,Neib,newT1,newQuad,NC,QFPossible)
		Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)
!Part 6:
        !------------------------------- Modify split Front Info -------------------------------
	    print*,'Smaller Front Removed!!! ',CF,' : ',FrontEdges(CF,LeftVertex),FrontEdges(CF,RightVertex)
	    States(CF) = Processed

	    if(FrontEdges(NFI,LeftVertex) == NK) then

		    FrontEdges(NFI,LeftVertex) = N1
		    Fronts = Fronts + 1
		    FrontEdges(Fronts,LeftVertex) = NV
		    FrontEdges(Fronts,RightVertex) = N1
		    FrontEdges(Fronts,Level) = FrontEdges(NFI,Level)
		    States(Fronts) = Recently_Added
		    print*,'Front Added!!! ',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)
		
	    else

		    FrontEdges(NFI,RightVertex) = N1
		    Fronts = Fronts + 1
		    FrontEdges(Fronts,LeftVertex) = N1
		    FrontEdges(Fronts,RightVertex) = NV
		    FrontEdges(Fronts,Level) = FrontEdges(NFI,Level)
		    States(Fronts) = Recently_Added
		    print*,'Front Added!!! ',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

        endif
        
        !---------------------------------- Updating NFI element -------------------------------
        
        FrontEdges(NFI,Element) = getNeibour(Dim,Corn,Neib,NF_Vertex,N1,newQ)
        
		!--------------------------- Specifying new Front's Element ---------------------------
		do I=1,4
			if(Corn(NC,I) == NV) I1 = I
			if(Corn(NC,I) == N1) I2 = I
		end do

		if(I1 < I2) then
			if(I2 - I1 == 1) then
				FrontEdges(Fronts,Element) = Neib(NC,I1)
			else
				FrontEdges(Fronts,Element) = Neib(NC,I2)
			endif
		else
			if(I1 - I2 == 1) then
				FrontEdges(Fronts,Element) = Neib(NC,I2)
			else
				FrontEdges(Fronts,Element) = Neib(NC,I1)
			endif
        endif
        
        !--------------------------------- Smoothing Q -----------------------------------------
	    newQuad(NK_i) = Corn(Q,NK_i)
	    newQuad(NF_Vertex_i) = Corn(Q,NF_Vertex_i)
	    newQuad(C_i) = Corn(Q,C_i)
	    newQuad(P_i) = Corn(Q,P_i)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Q,X,Y)

	    !------------------------------- Smoothing newQ ----------------------------------------
	    newQuad(NK_i) = Corn(newQ,NK_i)
	    newQuad(NF_Vertex_i) = Corn(newQ,NF_Vertex_i)
	    newQuad(C_i) = Corn(newQ,C_i)
	    newQuad(P_i) = Corn(newQ,P_i)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,newQ,X,Y)
    
    else !-------------- Case when edge recovery is not possible ----------
!Part 7:
        !-------------------------- Reconstructing Q element --------------------------------
        do I=1,4
            if(Corn(Q,I) == N2) then
                Corn(Q,I) = NF_Vertex
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,Q,C,NF_Vertex,E1)
        Call setNeibour(Dim,Corn,Neib,Q,NK,NF_Vertex,NFE)
        !-------------------------- Reconstructing NFE element ------------------------------
        do I=1,3
            if(Corn(NFE,I) == N1) then
                Corn(NFE,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,NFE,NK,NF_Vertex,Q)
        Call setNeibour(Dim,Corn,Neib,NFE,NK,NFEC,E2)
        
        !-------------------------- Modifying E1,E2 neibours -------------------------------
        
        Call setNeibour(Dim,Corn,Neib,E1,C,NF_Vertex,Q)
        Call setNeibour(Dim,Corn,Neib,E2,NK,NFEC,NFE)
        
        !-------------------------- Deleting newT1,newT2,newQ ------------------------------
        
        Call DeleteCell(Dim,newT1,Corn)
        Call DeleteCell(Dim,newT2,Corn)
        Call DeleteCell(Dim,newQ,Corn)
        
	endif

endif
!===========================================================================================
End Subroutine TransitionSplit
!*********************************************************************************************
