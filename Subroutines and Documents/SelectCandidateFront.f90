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
Subroutine SelectCandidateFront(Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,current_Level,Seeking_State,CF,of_State)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,Fronts,FrontEdges,States,X,Y,Seeking_State,current_Level
Intent(Out)::CF,of_State

Integer,Parameter::STATUS_THREE = 3
Integer,Parameter::Processed = -1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3

Integer::Dim,Fronts,Seeking_State,CF,of_State,C,I,J,K,L,M,index,current_Level
Integer,Dimension(1:Dim)::States,FList
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::Delayed
Real(8)::GetNorm,FN,N
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
CF = 0 !-------------------------------- Candidate Front -----------------------------------
C = 0  !----------------------- Counter for Fronts of the same level -----------------------

!Part 1:

do I=1,Fronts !---------------- Finding Fronts of the SAME LEVEL -----------------------
        
    if(States(I) /= Processed) then
            
        if(FrontEdges(I,Level) == current_Level) then !FrontEdges(I,Level) == L
        
            C = C + 1
            FList(C) = I    
            
        endif
            
    endif
        
end do

if(C /= 0) then
!Part 2:
    do I=1,C

		index = FList(I)
               
		if(States(index) == Seeking_State) then
			CF = index
			FN = GetNorm(Dim,FrontEdges(index,LeftVertex),FrontEdges(index,RightVertex),X,Y)
			of_State = States(CF)
			exit
        endif
        
    end do
!Part 3:	 
	do I=index,C

        if(States(I) == Seeking_State) then
            
			J = FList(I)

			N = GetNorm(Dim,FrontEdges(J,LeftVertex),FrontEdges(J,RightVertex),X,Y)

			if(FN > N) then
				
				FN = N
				CF = J 
				of_State = States(CF)
                
			endif
        
        endif

    end do

endif
!===========================================================================================
End Subroutine SelectCandidateFront
!*********************************************************************************************
