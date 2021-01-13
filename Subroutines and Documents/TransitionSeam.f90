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
Subroutine TransitionSeam(Dim,NC,NP,NBC,NK,NV,NF_Vertex,CF,NFI,ME,NFE,MFEC,NFEC,Fronts,Corn,Neib,FrontEdges,X,Y,States,BFP,MFN,NFN,TopEdgeRecoveryIsPossible)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBC,NK,NV,NF_Vertex,CF,NFI,ME,NFE,MFEC,NFEC,BFP,MFN,NFN
Intent(InOut)::Corn,Neib,FrontEdges,X,Y,States,Fronts,TopEdgeRecoveryIsPossible

Integer,Parameter::Recently_Added = -2
Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4
Integer,Parameter::PF=2	!---- PF(Processed Fronts)

Integer::Dim,NC,NP,NBC,NK,NV,NF_Vertex,CF,NFI,ME,NFE,MFEC,NFEC,Fronts,N1,I,I1,I2,J1,J2,J3,Q,newT1,newT2,C,C_i,P_i,NV_i,NF_Vertex_i,NK_i,E1,E2,getNeibour
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::TopEdgeRecoveryIsPossible,QFPossible
Real(8)::MFN,NFN
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!____________________________ Case when ME and its adjacent Quad are split _________________

if(MFN > NFN) then
!Part 1:    
	!-------------------------- Finding the Quad (the Splitting Quad) ---------------------- 
	do I=1,3
		if(Corn(ME,I) == MFEC) then
			Q = Neib(ME,I)
			exit
		endif
	end do

	!----------------------- Finding the other Adjacent Corner to NK in Q ------------------
	do I=1,4
		if(Corn(Q,I) == NK) then
			NK_i = I
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

	if(Corn(Q,I1) /= NV) then
		C = Corn(Q,I1)
		C_i = I1
		NV_i = I2
	else
		C = Corn(Q,I2)
		C_i = I2
		NV_i = I1
	endif
	
	!------------------- Specifying Index of the other Corner of splitting Quad ------------
	do I=1,4
		if(I/=NK_i .And. I/=NV_i .And. I/=C_i) then
			P_i = I
			exit
		endif
	end do
	!----------------------- Specifying Indices of ME corners ------------------------------
	do I=1,3
		if(Corn(ME,I) == MFEC) J1 = I !---------- Corresponding Index of MFEC --------------
		if(Corn(ME,I) == NV)   J2 = I !---------- Corresponding Index of NV ----------------
		if(Corn(ME,I) == NK)   J3 = I !---------- Corresponding Index of NK ----------------
    end do
!Part 2:    
	!--------------------------- N1: Mid Point of Larger Front -----------------------------
	NP = NP + 1
	N1 = NP
	X(N1) = (X(NK)+X(NV))/2
	Y(N1) = (Y(NK)+Y(NV))/2
	!--------------------- Composing new Triangle newT1(N1,NK,C) ---------------------------
	NC = NC + 1
	newT1 = NC
	Corn(newT1,J1) = C
	Corn(newT1,J2) = NK
	Corn(newT1,J3) = N1
	Corn(newT1,4) = 0
	!--------------------- Composing new Triangle newT2(N1,MFEC,NK) ------------------------
	NC = NC + 1
	newT2 = NC
	Corn(newT2,J1) = MFEC
	Corn(newT2,J2) = N1
	Corn(newT2,J3) = NK
	Corn(newT2,4) = 0
	!----------------------------- Determining newT1 neibours ------------------------------
    
    Call setNeibour(Dim,Corn,Neib,newT1,NK,N1,newT2)
    Call setNeibour(Dim,Corn,Neib,newT1,C,N1,Q)
    E1 = getNeibour(Dim,Corn,Neib,NK,C,Q)
    Call setNeibour(Dim,Corn,Neib,newT1,NK,C,E1)
	Neib(newT1,4) = 0
	
	!--------------------------- Updating Neibours of Neibours of newT1 --------------------
    
    Call setNeibour(Dim,Corn,Neib,E1,NK,C,newT1)
    
	!----------------------------- Determining newT2 neibours ------------------------------
    
    Call setNeibour(Dim,Corn,Neib,newT2,NK,N1,newT1)
    E2 = getNeibour(Dim,Corn,Neib,NK,MFEC,ME)
    Call setNeibour(Dim,Corn,Neib,newT2,NK,MFEC,E2)
    Call setNeibour(Dim,Corn,Neib,newT2,N1,MFEC,ME) 
	Neib(newT2,4) = 0
    
	!--------------------------- Updating Neibours of Neibours of newT2 --------------------
	
    Call setNeibour(Dim,Corn,Neib,E2,NK,MFEC,newT2)
  
	!--------------------------- Modifying the Splitting Quad ------------------------------
	do I=1,4
		if(Corn(Q,I) == NK) then 
            Corn(Q,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,Q,N1,C,newT1)
    
	!------------------------------- Modifying ME triangle ---------------------------------
	do I=1,3
		if(Corn(ME,I) == NK) then
            Corn(ME,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,ME,N1,MFEC,newT2)

	!-------------------- Composing the newQ(NF_vertex,N1,C,NK) Quad -----------------------

!Part 3:    
    
	Call TopEdgeRecovery(Dim,NC,NF_vertex,N1,Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
    
	if(TopEdgeRecoveryIsPossible) then
!Part 4:        
		newQuad(NK_i) = NK
		newQuad(NV_i) = NF_Vertex
		newQuad(C_i) = C
		newQuad(P_i) = N1
		Call QuadrilateralFormation(Dim,Corn,Neib,newT1,newQuad,NC,QFPossible)
		Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)

	    !------------------------------------ Smoothing Q --------------------------------------

	    newQuad(1) = Corn(Q,1)
	    newQuad(2) = Corn(Q,2)
	    newQuad(3) = Corn(Q,3)
	    newQuad(4) = Corn(Q,4)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Q,X,Y)
!Part 5:
	    !------------------------------- Modify split Front Info -------------------------------
	
	    if(FrontEdges(CF,RightVertex) == NK) then

		    FrontEdges(CF,RightVertex) = N1

			print *,'Smaller Front Removed!!! : ',NFI,' : ',FrontEdges(NFI,LeftVertex),FrontEdges(NFI,RightVertex)
			States(NFI) = Processed

			Fronts = Fronts + 1
			FrontEdges(Fronts,LeftVertex) = N1
			FrontEdges(Fronts,RightVertex) = NF_vertex
			FrontEdges(Fronts,Level) = FrontEdges(CF,Level)
			States(Fronts) = Recently_Added
			print *,'Front Added!!! : ',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

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
            
            !-------------------- Modifying CF Element --------------------- 
            
            do I=1,4
				if(Corn(Q,I) == NV) I1 = I
				if(Corn(Q,I) == N1) I2 = I
			end do

			if(I1 < I2) then
				if(I2 - I1 == 1) then
					FrontEdges(CF,Element) = Neib(Q,I1)
				else
					FrontEdges(CF,Element) = Neib(Q,I2)
				endif
			else
				if(I1 - I2 == 1) then
					FrontEdges(CF,Element) = Neib(Q,I2)
				else
					FrontEdges(CF,Element) = Neib(Q,I1)
				endif
            endif
			
	    else

		    FrontEdges(CF,LeftVertex) = N1

			print *,'Smaller Front Removed!!! : ',NFI,' : ',FrontEdges(NFI,LeftVertex),FrontEdges(NFI,RightVertex)
			States(NFI) = Processed

			Fronts = Fronts + 1
			FrontEdges(Fronts,LeftVertex) = NF_vertex
			FrontEdges(Fronts,RightVertex) = N1
			FrontEdges(Fronts,Level) = FrontEdges(CF,Level)
			States(Fronts) = Recently_Added
			print *,'Front Added!!! : ',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

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
            
            !-------------------- Modifying CF Element --------------------- 
            
            do I=1,4
				if(Corn(Q,I) == NV) I1 = I
				if(Corn(Q,I) == N1) I2 = I
			end do

			if(I1 < I2) then
				if(I2 - I1 == 1) then
					FrontEdges(CF,Element) = Neib(Q,I1)
				else
					FrontEdges(CF,Element) = Neib(Q,I2)
				endif
			else
				if(I1 - I2 == 1) then
					FrontEdges(CF,Element) = Neib(Q,I2)
				else
					FrontEdges(CF,Element) = Neib(Q,I1)
				endif
            endif		

        endif
        
    else !-------------------- Case when top edge recovery is not possible -----------------
!Part 6:        
        !------------------------ Reconstructing Q element --------------------------------- 
        
        do I=1,4
            if(Corn(Q,I) == N1) then
                Corn(Q,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,Q,NK,C,E1)
        
        !------------------------ Reconstructing ME element --------------------------------
        
        do I=1,3
            if(Corn(ME,I) == N1) then
                Corn(ME,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,ME,NK,MFEC,E2)
        
        !-------------------------- Modifying E1,E2 neibours -------------------------------
        
        Call setNeibour(Dim,Corn,Neib,E1,NK,C,Q)
        Call setNeibour(Dim,Corn,Neib,E2,NK,MFEC,ME)
        
        !--------------------------- Deleting newT1,newT2 ----------------------------------
        
        Call DeleteCell(Dim,newT1,Corn)
        Call DeleteCell(Dim,newT2,Corn)
        
    endif

!___________________________ Case when NFE and its adjacent Quad are split _________________	 		
else
!Part 1:    
	!-------------------------- Finding the Quad (the Splitting Quad) ---------------------- 
	do I=1,3
		if(Corn(NFE,I) == NFEC) then
			Q = Neib(NFE,I)
			exit
		endif
	end do

	!----------------------- Finding the other Adjacent Corner to NK in Q ------------------
	do I=1,4
		if(Corn(Q,I) == NK) then
			NK_i = I
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

	if(Corn(Q,I1) /= NF_vertex) then
		C = Corn(Q,I1)
		C_i = I1
		NF_vertex_i = I2
	else
		C = Corn(Q,I2)
		C_i = I2
		NF_vertex_i = I1
	endif

	!------------------- Specifying Index of the other Corner of splitting Quad ------------
	do I=1,4
		if(I/=NK_i .And. I/=NF_vertex_i .And. I/=C_i) then
			P_i = I
			exit
		endif
	end do
	!----------------------- Specifying Indices of NFE corners ------------------------------
	do I=1,3
		if(Corn(NFE,I) == NFEC)      J1 = I !-------- Corresponding Index of NFEC -----------
		if(Corn(NFE,I) == NF_vertex) J2 = I !-------- Corresponding Index of NF_vertex ------
		if(Corn(NFE,I) == NK)        J3 = I !-------- Corresponding Index of NK -------------
    end do
!Part 2:    
	!--------------------------- N1: Mid Point of Larger Front -----------------------------
	NP = NP + 1
	N1 = NP
	X(N1) = (X(NK)+X(NF_vertex))/2
	Y(N1) = (Y(NK)+Y(NF_vertex))/2
	!--------------------- Composing new Triangle newT1(N1,NK,C) ---------------------------
	NC = NC + 1
	newT1 = NC
	Corn(newT1,J1) = C
	Corn(newT1,J2) = NK
	Corn(newT1,J3) = N1
	Corn(newT1,4) = 0
	!--------------------- Composing new Triangle newT2(N1,MFEC,NK) ------------------------
	NC = NC + 1
	newT2 = NC
	Corn(newT2,J1) = NFEC
	Corn(newT2,J2) = N1
	Corn(newT2,J3) = NK
	Corn(newT2,4) = 0
	!----------------------------- Determining newT1 neibours ------------------------------
    
    Call setNeibour(Dim,Corn,Neib,newT1,NK,N1,newT2)
    Call setNeibour(Dim,Corn,Neib,newT1,C,N1,Q)
    E1 = getNeibour(Dim,Corn,Neib,NK,C,Q)
    Call setNeibour(Dim,Corn,Neib,newT1,C,NK,E1)
	Neib(newT1,4) = 0
	
	!--------------------------- Updating Neibours of Neibours of newT1 --------------------
    
    Call setNeibour(Dim,Corn,Neib,E1,C,NK,newT1)
    
	!----------------------------- Determining newT2 neibours ------------------------------
	Call setNeibour(Dim,Corn,Neib,newT2,NK,N1,newT1)
    E2 = getNeibour(Dim,Corn,Neib,NK,NFEC,NFE)
    Call setNeibour(Dim,Corn,Neib,newT2,NK,NFEC,E2)
    Call setNeibour(Dim,Corn,Neib,newT2,N1,NFEC,NFE)
	Neib(newT2,4) = 0
	!--------------------------- Updating Neibours of Neibours of newT2 --------------------
    
    Call setNeibour(Dim,Corn,Neib,E2,NK,NFEC,newT2)

	!--------------------------- Modifying the Splitting Quad ------------------------------
	do I=1,4
		if(Corn(Q,I) == NK) then
            Corn(Q,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,Q,N1,C,newT1)
	!------------------------------- Modifying NFE triangle ---------------------------------
	do I=1,3
		if(Corn(NFE,I) == NK) then
            Corn(NFE,I) = N1
            exit
        endif
    end do
    
    Call setNeibour(Dim,Corn,Neib,NFE,N1,NFEC,newT2)
	!-------------------- Composing the newQ(NF_vertex,N1,C,NK) Quad -----------------------
!Part 3: 
	Call TopEdgeRecovery(Dim,NC,NV,N1,Corn,Neib,X,Y,FrontEdges,Fronts,States,TopEdgeRecoveryIsPossible)
    
	if(TopEdgeRecoveryIsPossible) then
!Part 4:         
		newQuad(NK_i) = NK
		newQuad(NF_vertex_i) = NV
		newQuad(C_i) = C
		newQuad(P_i) = N1
		Call QuadrilateralFormation(Dim,Corn,Neib,newT1,newQuad,NC,QFPossible)
		Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,NC,X,Y)

	    !------------------------------------ Smoothing Q --------------------------------------

	    newQuad(1) = Corn(Q,1)
	    newQuad(2) = Corn(Q,2)
	    newQuad(3) = Corn(Q,3)
	    newQuad(4) = Corn(Q,4)
	    Call Smooth(Dim,NC,NBC,BFP,Fronts,Corn,Neib,FrontEdges,States,newQuad,Q,X,Y)
!Part 5:
	    !------------------------------- Modify split Front Info -------------------------------

	    if(FrontEdges(NFI,RightVertex) == NK) then
		
		    FrontEdges(NFI,RightVertex) = N1
			
			print *,'Smaller Front Removed!!! :',CF,' : ',FrontEdges(CF,LeftVertex),FrontEdges(CF,RightVertex)
			States(CF) = Processed

			Fronts = Fronts + 1
			FrontEdges(Fronts,LeftVertex) = N1
			FrontEdges(Fronts,RightVertex) = NV
			FrontEdges(Fronts,Level) = FrontEdges(NFI,Level)
			States(Fronts) = Recently_Added
			print *,'Front Added!!! :',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

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
            
            !-------------------- Modifying NFI Element --------------------- 
            
            do I=1,4
				if(Corn(Q,I) == NF_vertex) I1 = I
				if(Corn(Q,I) == N1) I2 = I
			end do

			if(I1 < I2) then
				if(I2 - I1 == 1) then
					FrontEdges(NFI,Element) = Neib(Q,I1)
				else
					FrontEdges(NFI,Element) = Neib(Q,I2)
				endif
			else
				if(I1 - I2 == 1) then
					FrontEdges(NFI,Element) = Neib(Q,I2)
				else
					FrontEdges(NFI,Element) = Neib(Q,I1)
				endif
            endif

	    else

		    FrontEdges(NFI,LeftVertex) = N1
			
			print *,'Smaller Front Removed!!! :',CF,' : ',FrontEdges(CF,LeftVertex),FrontEdges(CF,RightVertex)
			States(CF) = Processed

			Fronts = Fronts + 1
			FrontEdges(Fronts,LeftVertex) = NV
			FrontEdges(Fronts,RightVertex) = N1
			FrontEdges(Fronts,Level) = FrontEdges(NFI,Level)
			States(Fronts) = Recently_Added
			print *,'Front Added!!! :',Fronts,' : ',FrontEdges(Fronts,LeftVertex),FrontEdges(Fronts,RightVertex)

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
            
            !-------------------- Modifying NFI Element --------------------- 
            
            do I=1,4
				if(Corn(Q,I) == NF_vertex) I1 = I
				if(Corn(Q,I) == N1) I2 = I
			end do

			if(I1 < I2) then
				if(I2 - I1 == 1) then
					FrontEdges(NFI,Element) = Neib(Q,I1)
				else
					FrontEdges(NFI,Element) = Neib(Q,I2)
				endif
			else
				if(I1 - I2 == 1) then
					FrontEdges(NFI,Element) = Neib(Q,I2)
				else
					FrontEdges(NFI,Element) = Neib(Q,I1)
				endif
            endif

        endif
        
    else !---------------------- Case when edge recovery is impossible ---------------------
!Part 6:         
        !------------------------ Reconstructing Q element --------------------------------- 
        
        do I=1,4
            if(Corn(Q,I) == N1) then
                Corn(Q,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,Q,NK,C,E1)
        
        !------------------------ Reconstructing ME element --------------------------------
        
        do I=1,3
            if(Corn(NFE,I) == N1) then
                Corn(NFE,I) = NK
                exit
            endif
        end do
        
        Call setNeibour(Dim,Corn,Neib,NFE,NK,NFEC,E2)
        
        !-------------------------- Modifying E1,E2 neibours -------------------------------
        
        Call setNeibour(Dim,Corn,Neib,E1,NK,C,Q)
        Call setNeibour(Dim,Corn,Neib,E2,NK,NFEC,NFE)
        
        !--------------------------- Deleting newT1,newT2 ----------------------------------
        
        Call DeleteCell(Dim,newT1,Corn)
        Call DeleteCell(Dim,newT2,Corn)
        
    endif

endif
!===========================================================================================
End Subroutine TransitionSeam
!*********************************************************************************************
