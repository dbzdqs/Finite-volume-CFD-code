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
Subroutine FiveThreeOppoFive(Dim,NC,NP,Corn,Neib,X,Y,NBE,BFP,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::NC,NP,Corn,Neib,X,Y,done

Integer::Dim,NC,NP,NBE,TEC,ME,I,J,A,Ai,B,Bi,C,Ci,D,Di,Na,Nb,Nc1,Nd,P1,Np1,Elm,getNeibour,Nab
Integer,Dimension(1:1000)::QElmsA,TElms,QElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,areAdjacent
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do I=1,4

	done = .False.

	Ai = I
	A = Corn(ME,Ai)

	Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElmsA,TEC,Na)
!Part 1:
	if(Na == 3) then !------------------------ Condition on A Valance ----------------------

		Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

		Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc1)

		if(Nc1 == 5) then !--------------------- Condition on C Valance --------------------

			!---------------------- Finding Next and Previous Corner to 'A' ----------------

			Call getNextCorner(Dim,Corn,ME,A,B,Bi)
			Call getPrevCorner(Dim,Corn,ME,A,D,Di)

			Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElms,TEC,Nb)
			Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)

			if(Nb == 4 .And. Nd == 4) then !------- Checking 'B' and 'D' Valances ----------

				!------------------ Finding Nab Element (ME naibours) ----------------------

				Nab = getNeibour(Dim,Corn,Neib,A,B,ME)

				do J=1,4

					if(areAdjacent(Dim,Corn,A,Corn(Nab,J),Nab) .And. Corn(Nab,J) /= B) then

						P1 = Corn(Nab,J)
						exit

					endif

				end do

				Call GetSurroundingElements(Dim,Corn,Neib,P1,Nab,TElms,QElms,TEC,Np1)

				if(Np1 == 5) then
!Part 2:
					!----- Collapse ME

					Call CollapseOperation(Dim,Corn,Neib,X,Y,NBE,BFP,D,A,B,ME,done)
					Call NodeElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,Nab,done)

					if(done) exit

				endif

			endif !----------------- End of Checking 'B' and 'D' Valances ------------------

		endif !------------------------ End of Condition on C Valance ----------------------

	endif !------------------------- End of Condition on A Valance -------------------------

end do !-------------------------- End of Main Loop ('I' Loop) -----------------------------
!===========================================================================================
End Subroutine FiveThreeOppoFive
!*********************************************************************************************
