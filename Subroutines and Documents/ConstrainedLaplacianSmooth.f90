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
Subroutine ConstrainedLaplacianSmooth(Dim,NBC,BEP,Corn,Neib,X,Y,P,E)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBC,BEP,Corn,Neib,P,E
Intent(InOut)::X,Y

Integer::Dim,NBC,P,E,TEC,QEC,I,J,K,Count,Q1,Q2
Integer,Dimension(1:1000)::TElms,QElms,PList,TEF,QEF
Integer,Dimension(1:Dim,1:2)::BEP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent,IsInTheList,PolygonIsConvex,IsOnExteriorBoundary,ElementInverted
Real(8)::earlyX,earlyY,Cx,Cy,QuadQuality
Real(8),Dimension(1:100)::PreBetaValues,PostBetaValues
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

if(.Not. IsOnExteriorBoundary(Dim,NBC,BEP,P)) then

!Part 1:
    
	Call GetSurroundingElements(Dim,Corn,Neib,P,E,TElms,QElms,TEC,QEC)
	Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms,QEC,TElms,TEC,TEF,QEF)

!Part 2:
    
	if(QEC > 2) then !----------->>> Notice: Since a polygon has 3 or more vertices ----------------

		Count = 0

		earlyX = X(P)
		earlyY = Y(P)

		do I=1,QEC !------------- Finding Corners connected to P Consecutively ---------------------

			Q1 = QElms(I)

			if(I+1 > QEC) then
				
				Q2 = QElms(1)

			else

				Q2 = QElms(I+1)

			endif

			do J=1,4 !---- Finding other point of the common edge between Q1 and Q2 that has P -----

				do K=1,4

					if(Corn(Q1,J) == Corn(Q2,K) .And. Corn(Q1,J) /= P .And. areAdjacent(Dim,Corn,Corn(Q1,J),P,Q1)) then
						if(.Not. IsInTheList(PList,Count,Corn(Q1,J))) then
							Count = Count + 1
							PList(Count) = Corn(Q1,J)
						endif
					endif

				end do

			end do

		end do
!Part 3:
		do I=1,QEC

			PreBetaValues(I) = QuadQuality(Dim,X,Y,Corn(QElms(I),1),Corn(QElms(I),2),Corn(QElms(I),3),Corn(QElms(I),4))

        end do

!Part 4:
        
		Call CenterOfPolygon(Dim,Count,PList,X,Y,Cx,Cy)
			
		X(P) = Cx
		Y(P) = Cy
		
!Part 5:
        
		do I=1,QEC

			PostBetaValues(I) = QuadQuality(Dim,X,Y,Corn(QElms(I),1),Corn(QElms(I),2),Corn(QElms(I),3),Corn(QElms(I),4))

        end do
	
!Part 6:        
        
		if(ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)) then
			
			X(P) = earlyX
			Y(P) = earlyY
			
		else

			do I=1,QEC
				
				if(PostBetaValues(I) < PreBetaValues(I)) then

					X(P) = earlyX
					Y(P) = earlyY
					exit

				endif

			end do
			
		endif 

	endif

endif

!===========================================================================================
End Subroutine ConstrainedLaplacianSmooth 
!*********************************************************************************************
