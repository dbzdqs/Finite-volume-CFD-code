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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine MergeTriangles(Dim,NC,NP,Corn,Neib,X,Y,FrontEdges,States,FrontIndex,Fronts,RV,LV,newQuad)
Implicit None
!===========================================================================================
Intent(In)::Dim,States,FrontIndex,Fronts,RV,LV
Intent(InOut)::X,Y,Corn,Neib,FrontEdges,NP,NC,newQuad

Integer,Parameter::to_Left = 1
Integer,Parameter::to_Right = 2
Integer,Parameter::Element=4

Integer::Dim,NC,NP,RV,LV,MFEC,ME,I,NR,NL,NML,NMR,Nn,Fronts,FrontIndex
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::isOnTheBoundary,NumOfEdgesInTheLoopIsEven
Real(8)::GetAngle
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
print *,'>>>>>>>>>>>>>>>>> Merging Triangles <<<<<<<<<<<<<<<<<<<'

!Part 1:

ME = FrontEdges(FrontIndex,Element)

do I=1,3
	if(Corn(ME,I) == RV) then
		NL = Neib(ME,I)
	endif
	if(Corn(ME,I) == LV) then
		NR = Neib(ME,I)
	endif
	if(Corn(ME,I)/=RV .And. Corn(ME,I)/=LV) then
		MFEC = Corn(ME,I)
	endif
end do

if(isOnTheBoundary(Dim,Fronts,MFEC,FrontEdges,States)) then
!Part 2:	
	if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,RV,MFEC,to_Left)) then

		do I=1,3
			if(Neib(NL,I) == ME) then
				NML = Corn(NL,I)
				exit
			endif
		end do

		if(isOnTheBoundary(Dim,Fronts,NML,FrontEdges,States)) then 

			if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LV,NML,to_Left)) then
				
				newQuad(1) = MFEC
				newQuad(2) = NML

			else

				Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NL,Nn)
				newQuad(1) = MFEC
				newQuad(2) = Nn

			endif

		else

			newQuad(1) = MFEC
			newQuad(2) = NML

		endif

	else
!Part 3:
		do I=1,3
			if(Neib(NR,I) == ME) then
				NMR = Corn(NR,I)
				exit
			endif
		end do

		if(isOnTheBoundary(Dim,Fronts,NMR,FrontEdges,States)) then

			if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,RV,NMR,to_Right)) then

				newQuad(1) = NMR
				newQuad(2) = MFEC

			else

				Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NR,Nn)
				newQuad(1) = Nn
				newQuad(2) = MFEC

			endif

		else

			newQuad(1) = NMR
			newQuad(2) = MFEC

		endif

	endif

else
!Part 4:	
	if(Corn(NL,4) == 0 .And. Corn(NR,4) == 0) then
		
		do I=1,3
			if(Neib(NL,I) == ME) then
				NML = Corn(NL,I)
			endif
			if(Neib(NR,I) == ME) then
				NMR = Corn(NR,I)
			endif
		end do

		if(GetAngle(Dim,RV,LV,NMR,X,Y) < GetAngle(Dim,LV,RV,NML,X,Y)) then
			
			if(isOnTheBoundary(Dim,Fronts,NMR,FrontEdges,States)) then
				
				if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,RV,NMR,to_Right)) then
					
					newQuad(1) = NMR
					newQuad(2) = MFEC

				else

					if(isOnTheBoundary(Dim,Fronts,NML,FrontEdges,States)) then

						if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LV,NML,to_Left)) then
							
							newQuad(1) = MFEC
							newQuad(2) = NML

						else

							Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NR,Nn)
							newQuad(1) = Nn
							newQuad(2) = MFEC

						endif

					else

						newQuad(1) = MFEC
						newQuad(2) = NML

					endif

				endif

			else

				newQuad(1) = NMR
				newQuad(2) = MFEC

			endif

		else

			if(isOnTheBoundary(Dim,Fronts,NML,FrontEdges,States)) then
				
				if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LV,NML,to_Left)) then
					
					newQuad(1) = MFEC
					newQuad(2) = NML

				else

					if(isOnTheBoundary(Dim,Fronts,NMR,FrontEdges,States)) then
						
						if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,RV,NMR,to_Right)) then
							
							newQuad(1) = NMR
							newQuad(2) = MFEC

						else

							Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NL,Nn)
							newQuad(1) = MFEC
							newQuad(2) = Nn

						endif

					else

						newQuad(1) = NMR
						newQuad(2) = MFEC

					endif

				endif

			else

				newQuad(1) = MFEC
				newQuad(2) = NML

			endif

		endif
  
	elseif(Corn(NL,4) == 0 .And. Corn(NR,4) /= 0) then
!Part 5:
		do I=1,3
			if(Neib(NL,I) == ME) then
				NML = Corn(NL,I)
				exit
			endif
		end do

		if(isOnTheBoundary(Dim,Fronts,NML,FrontEdges,States)) then
			
			if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,LV,NML,to_Left)) then	
				
				newQuad(1) = MFEC
				newQuad(2) = NML

			else

				Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NL,Nn)
				newQuad(1) = MFEC
				newQuad(2) = Nn

			endif

		else

			newQuad(1) = MFEC
			newQuad(2) = NML

		endif

	elseif(Corn(NL,4) /= 0 .And. Corn(NR,4) == 0) then
!Part 6:		
		do I=1,3
			if(Neib(NR,I) == ME) then
				NMR = Corn(NR,I)
				exit
			endif
		end do

		if(isOnTheBoundary(Dim,Fronts,NMR,FrontEdges,States)) then
			
			if(NumOfEdgesInTheLoopIsEven(Dim,Corn,Neib,FrontEdges,Fronts,States,FrontIndex,RV,NMR,to_Right)) then	
				
				newQuad(1) = NMR
				newQuad(2) = MFEC

			else

				Call TriangleSplit(Dim,Corn,Neib,X,Y,States,NC,NP,FrontEdges,Fronts,NR,Nn)
				newQuad(1) = Nn
				newQuad(2) = MFEC

			endif

		else

			newQuad(1) = NMR
			newQuad(2) = MFEC

		endif

	endif
	
endif
!===========================================================================================
End Subroutine MergeTriangles
!*********************************************************************************************
