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
Subroutine TwoCollapsesOperation(Dim,Corn,Neib,X,Y,NBE,BFP,ME,done)	! fig. 8
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::Corn,Neib,X,Y,done

Integer::Dim,NBE,TEC,ME,NE,I,J,A,Ai,B,Bi,C,Ci,D,Di,Na,Nb,Nc1,Nd,P0,P1,P2,P3,Pi,Np2,Np3,getNeibour
Integer,Dimension(1:1000)::QElmsA,TElms,QElms
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,areAdjacent,ElementFound
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!--->>> Reference: Canann, S. A., S. Muthukrishnan and B. Phillips (1994),"Topological improvement procedures for quadrilateral and triangular finite element meshes"

do I=1,4

	done = .False.
	ElementFound = .False.

	Ai = I
	A = Corn(ME,Ai)
!Part 1:
	Call GetSurroundingElements(Dim,Corn,Neib,A,ME,TElms,QElmsA,TEC,Na)

	if(Na == 5) then !------------------------ Condition on A Valance ----------------------

		Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

		Call GetSurroundingElements(Dim,Corn,Neib,C,ME,TElms,QElms,TEC,Nc1)

		if(Nc1 == 4) then !--------------------- Condition on C Valance --------------------

			Call getNextCorner(Dim,Corn,ME,A,B,Bi)
			Call getPrevCorner(Dim,Corn,ME,A,D,Di)

			Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElms,TEC,Nb)
			Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElms,TEC,Nd)

			if(Nb == 3 .And. Nd == 4) then !------- Checking 'B' and 'D' Valances ----------
			
				P0 = B
				P1 = D
				ElementFound = .True.

			elseif(Nb == 4 .And. Nd == 3) then

				P0 = D
				P1 = B
				ElementFound = .True.

			endif !----------------- End of Checking 'B' and 'D' Valances ------------------

			if(ElementFound) then

				NE = getNeibour(Dim,Corn,Neib,C,P1,ME)

				do J=1,4

					if(areAdjacent(Dim,Corn,P1,Corn(NE,J),NE) .And. Corn(NE,J) /= C) then

						P2 = Corn(NE,J)
						exit

					endif

				end do

				Call GetSurroundingElements(Dim,Corn,Neib,P2,NE,TElms,QElms,TEC,Np2)

				Call getOppoCorner(Dim,Corn,NE,P1,P3,Pi)
				Call GetSurroundingElements(Dim,Corn,Neib,P3,NE,TElms,QElms,TEC,Np3)

				if(Np2 == 3 .And. Np3 == 5) then 
!Part 2:
					Call CollapseOperation(Dim,Corn,Neib,X,Y,NBE,BFP,P2,P3,C,NE,done)
					Call CollapseOperation(Dim,Corn,Neib,X,Y,NBE,BFP,P0,A,P1,ME,done)

					if(done) exit

				endif

			endif

		endif !------------------------ End of Condition on C Valance ----------------------

	endif !------------------------- End of Condition on A Valance -------------------------

end do !-------------------------- End of Main Loop ('I' Loop) -----------------------------
!===========================================================================================
End Subroutine TwoCollapsesOperation
!*********************************************************************************************
