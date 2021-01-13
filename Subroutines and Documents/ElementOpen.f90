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
Subroutine ElementOpen(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,ME,NBE,BFP
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,ME,E1,E2,I,J,K,TEC,Elm,A,Ai,A2,B,C,Ci,Na,P1,Pi1,P2,Pi2,P3,Np1,Np2,index,getNeibour,Cnt
Integer,Dimension(1:1000)::QElmsA,TElms,QElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::Element_Found,done,areAdjacent
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
do I=1,4

    done = .False.
	Element_Found = .False.

	Ai = I
	A = Corn(ME,Ai)

	Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElmsA,TEC,Na)

	if(Na >= 6) then !------------------------ Condition on A Valance ----------------------
!Part 1:
        !---------------- Finding a Corner (C) having valence of three or less -------------
        C = 0
        Call getNextCorner(Dim,Corn,ME,A,P1,Pi1)
        P2 = P1
        P3 = 0
		Elm = ME

		do
            Call GetSurroundingElements(Dim,Corn,Neib,P1,Elm,TElms,QElms,TEC,Np1)
            
            if(Np1 <= 3) then
                
                C = P1
                E1 = Elm
                exit
                
            elseif(Np1 == 4) then
            
                P3 = P1 !---- Alternative point in case when C couldn't be found ----
                E1 = Elm
                
            endif
            
			do J=1,4

				if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P1) then
							
					P1 = Corn(Elm,J)
                    Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)
					exit

				endif

			end do

			if(P1 == P2) exit

        end do
        
        if(C == 0) then !-------------------- In case C not found ----------------
            
            if(P3 /= 0) then
                C = P3    
            else
                C = P2
                E1 = ME
            endif
            
        endif
!Part 2:
		!-------------------------- Finding oposite corner to C -------------------------

		P1 = C
		Elm = E1
		Cnt = 0

		do

			do J=1,4

				if(areAdjacent(Dim,Corn,A,Corn(Elm,J),Elm) .And. Corn(Elm,J) /= P1) then

					Cnt = Cnt + 1
							
					P1 = Corn(Elm,J)

					if(Cnt /= Na/2) then

						Elm = getNeibour(Dim,Corn,Neib,A,P1,Elm)

					endif

					exit

				endif

			end do

			if(Cnt == Na/2) exit

        end do
        
        E2 = Elm
!Part 3:			
		Call ElementOpenningOperation(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,C,A,P1,E1,E2,done)

        if(done) then
		    exit    
        endif

	endif !----------------------- End of Condition on A Valance ---------------------------
														
end do !---------------------------- End of Main Loop ('I' Loop) ---------------------------

!===========================================================================================
End Subroutine ElementOpen
!*********************************************************************************************
