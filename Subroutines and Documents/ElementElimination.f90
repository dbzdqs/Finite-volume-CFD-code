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
Subroutine ElementElimination(Dim,NC,NBE,BFP,Corn,Neib,X,Y,E,ElementFound)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,NBE,BFP
Intent(InOut)::Corn,Neib,X,Y,ElementFound,E

Integer::Dim,NC,NBE,E,I,J,K,L,C1,C2,P,Q,Pi,Qi,NL1,NL2,NR1,NR2,TEC1,QEC1,TEC2,QEC2,I1,I2,PRIORITY
Integer,Dimension(1:1000)::TElmsC1,QElmsC1,TElmsC2,QElmsC2,TEF_C1,QEF_C1,TEF_C2,QEF_C2
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsOnExteriorBoundary,IsBoundaryElement,ElementInverted,IsBoundaryC1,IsBoundaryC2,ElementFound
Real(8)::x_value,y_value,x_C1,y_C1,x_C2,y_C2
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

ElementFound = .False.

if(.Not. IsBoundaryElement(Dim,Corn,NBE,BFP,E)) then
!Part 1:
	do J=1,2

		C1 = Corn(E,J)

		Call GetSurroundingElements(Dim,Corn,Neib,C1,E,TElmsC1,QElmsC1,TEC1,QEC1)

		if(QEC1 == 3 .And. .Not. IsOnExteriorBoundary(Dim,NBE,BFP,C1)) then

			K = J + 2
			C2 = Corn(E,K)

			Call GetSurroundingElements(Dim,Corn,Neib,C2,E,TElmsC2,QElmsC2,TEC2,QEC2)

			if(QEC2 == 3 .And. .Not. IsOnExteriorBoundary(Dim,NBE,BFP,C2)) then
!Part 2:
				Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElmsC1,QEC1,TElmsC1,TEC1,TEF_C1,QEF_C1)
				Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElmsC2,QEC2,TElmsC2,TEC2,TEF_C2,QEF_C2)

				x_C1 = X(C1)
				y_C1 = Y(C1)

				x_C2 = X(C2)
				y_C2 = Y(C2)
							
				x_value = (X(C1) + X(C2))/2
				y_value = (Y(C1) + Y(C2))/2

				X(C1) = x_value
				Y(C1) = y_value
				X(C2) = x_value
				Y(C2) = y_value
!Part 3:
				if(ElementInverted(Dim,Corn,X,Y,TElmsC1,QElmsC1,TEC1,QEC1,TEF_C1,QEF_C1) .Or. ElementInverted(Dim,Corn,X,Y,TElmsC2,QElmsC2,TEC2,QEC2,TEF_C2,QEF_C2)) then
							
					X(C1) = x_C1
					Y(C1) = y_C1

					X(C2) = x_C2
					Y(C2) = y_C2

				else
!Part 4:
					ElementFound = .True.

					if(J == 1) then

						Pi = J + 1
						Qi = K + 1

					elseif(J == 2) then
								
						Qi = J - 1
						Pi = K - 1

					endif

					P = Corn(E,Pi)
					Q = Corn(E,Qi)

					if(J > Pi) then
						if(J - Pi == 1) then
							NL1 = Neib(E, Pi)
						else
							NL1 = Neib(E, J)
						endif
					else
						if(Pi - J == 1) then
							NL1 = Neib(E, J)
						else
							NL1 = Neib(E, Pi)
						endif
					endif

					if(K > Pi) then
						if(K - Pi == 1) then
							NL2 = Neib(E, Pi)
						else
							NL2 = Neib(E, K)
						endif
					else
						if(Pi - K == 1) then
							NL2 = Neib(E, K)
						else
							NL2 = Neib(E, Pi)
						endif
					endif

					if(J > Qi) then
						if(J - Qi == 1) then
							NR1 = Neib(E, Qi)
						else
							NR1 = Neib(E, J)
						endif
					else
						if(Qi - J == 1) then
							NR1 = Neib(E, J)
						else
							NR1 = Neib(E, Qi)
						endif
					endif

					if(K > Qi) then
						if(K - Qi == 1) then
							NR2 = Neib(E, Qi)
						else
							NR2 = Neib(E, K)
						endif
					else
						if(Qi - K == 1) then
							NR2 = Neib(E, K)
						else
							NR2 = Neib(E, Qi)
						endif
					endif

					!========================================
!Part 5:
					do I=1,4
						if(Corn(NL1,I) == C1) I1 = I
						if(Corn(NL1,I) == P) I2 = I
					end do

					if(I1 > I2) then
						if(I1 - I2 == 1) then
							Neib(NL1,I2) = NL2
						else
							Neib(NL1,I1) = NL2
						endif
					else
						if(I2 - I1 == 1) then
							Neib(NL1,I1) = NL2
						else
							Neib(NL1,I2) = NL2
						endif
					endif

					do I=1,4
						if(Corn(NL2,I) == P) I1 = I
						if(Corn(NL2,I) == C2) I2 = I
					end do

					if(I1 > I2) then
						if(I1 - I2 == 1) then
							Neib(NL2,I2) = NL1
						else
							Neib(NL2,I1) = NL1
						endif
					else
						if(I2 - I1 == 1) then
							Neib(NL2,I1) = NL1
						else
							Neib(NL2,I2) = NL1
						endif
					endif

					do I=1,4
						if(Corn(NR1,I) == C1) I1 = I
						if(Corn(NR1,I) == Q) I2 = I
					end do

					if(I1 > I2) then
						if(I1 - I2 == 1) then
							Neib(NR1,I2) = NR2
						else
							Neib(NR1,I1) = NR2
						endif
					else
						if(I2 - I1 == 1) then
							Neib(NR1,I1) = NR2
						else
							Neib(NR1,I2) = NR2
						endif
					endif

					do I=1,4
						if(Corn(NR2,I) == C2) I1 = I
						if(Corn(NR2,I) == Q) I2 = I
					end do

					if(I1 > I2) then
						if(I1 - I2 == 1) then
							Neib(NR2,I2) = NR1
						else
							Neib(NR2,I1) = NR1
						endif
					else
						if(I2 - I1 == 1) then
							Neib(NR2,I1) = NR1
						else
							Neib(NR2,I2) = NR1
						endif
					endif
						

					!====================================================
!Part 6:						
					do I=1,QEC1
						do L=1,4
							if(Corn(QElmsC1(I),L) == C1) then
								Corn(QElmsC1(I),L) = C2	
							endif
						end do
					end do

					
					do I=1,4
						Corn(E,I) = Corn(E,1)
					end do
!Part 7:
					!---------->>>>>>> Smoothing neibours >>>>>>>>>>>>>>>>
					do I=1,4

						Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NL1,I),NL1)

					end do

					do I=1,4

						Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NL2,I),NL2)

					end do

					do I=1,4

						Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NR1,I),NR1)

					end do

					do I=1,4

						Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(NR2,I),NR2)

					end do

				endif

			endif

		endif
				
		if(ElementFound) exit

	end do

endif

!===========================================================================================
End Subroutine ElementElimination
!*********************************************************************************************
