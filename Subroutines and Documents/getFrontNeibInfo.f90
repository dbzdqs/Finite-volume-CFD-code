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
Subroutine getFrontNeibInfo(Dim,FrontEdges,Fronts,NK,NV,FrontIndex,Situation,Corn,Neib,States,NFE,NFI,NFEC,NF_Vertex)
Implicit None
!===========================================================================================
Intent(In)::Situation,Fronts,FrontEdges,Corn,Neib,States,NK,NV,FrontIndex
Intent(Out)::NFE,NFI,NFEC,NF_Vertex

Integer,Parameter::Processed=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Level=3
Integer,Parameter::Element=4

Integer::Dim,Fronts,NK,NV,Situation,I,J,ME,NFE,NFI,NFEC,NF_Vertex,FrontIndex,C,Prev,Next
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
!===========================================================================================
ME = FrontEdges(FrontIndex,Element)

Select Case(Situation)
							   
	Case(1) !---------------------------- In this case left Front Neibour information is gathered --------------------------------
!Part 1:	
		Prev = ME
		C = NV

		do
			do J=1,3
				if(Corn(Prev,J) == C) then
					Next = Neib(Prev,J)
					exit
				endif
			end do
			
			if(Next == 0) then
                exit    
            elseif(Corn(Next,4) /= 0) then
                exit    
            endif
			
			do J=1,3
				if(Corn(Prev,J)/=NK .And. Corn(Prev,J)/=C) then
					C = Corn(Prev,J)
					exit	
				endif
			end do

			Prev = Next
					
		end do
!Part 2:		
		NFE = Prev                      ! Neibouring Front Element(NFE)
		NFEC = C

		do J=1,3
			if(Corn(NFE,J)/= NK .And. Corn(NFE,J)/=NFEC) then
				NF_Vertex = Corn(NFE,J) ! Neibouring Front Vertex(NF_Vertex)
				exit	
			endif
		end do
!Part 3:		
		do I=1,Fronts
			if(States(I)/=Processed) then
				if(FrontEdges(I,LeftVertex) == NF_Vertex .And. FrontEdges(I,RightVertex) == NK) then           
					NFI = I			   ! Neibouring Front Index(NFI)
                    exit
				endif
			endif
		end do	
		
								            
	Case(2) ! ----------------------------- In this case right Front Neibour information is gathered -----------------------------
!Part 1:		
		Prev = ME
		C = NV

		do
			do J=1,3
				if(Corn(Prev,J) == C) then
					Next = Neib(Prev,J)
					exit
				endif
			end do
			
			if(Next == 0) then
                exit    
            elseif(Corn(Next,4) /= 0) then
                exit    
            endif
			
			do J=1,3
				if(Corn(Prev,J)/=NK .And. Corn(Prev,J)/=C) then
					C = Corn(Prev,J)
					exit	
				endif
			end do

			Prev = Next
					
		end do
!Part 2:		
		NFE = Prev                      ! Neibouring Front Element(NFE)
		NFEC = C

		do J=1,3
			if(Corn(NFE,J)/= NK .And. Corn(NFE,J)/=NFEC) then
				NF_Vertex = Corn(NFE,J) ! Neibouring Front Vertex(NF_Vertex)
				exit	
			endif
		end do
!Part 3:		
		do I=1,Fronts
			if(States(I)/=Processed) then
				if(FrontEdges(I,LeftVertex) == NK .And. FrontEdges(I,RightVertex) == NF_Vertex) then           
					NFI = I			   ! Neibouring Front Index(NFI)
                    exit
				endif
			endif
		end do

End Select
!===========================================================================================
End Subroutine getFrontNeibInfo
!*********************************************************************************************
