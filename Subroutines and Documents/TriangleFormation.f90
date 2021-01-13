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
Subroutine TriangleFormation(Dim,NC,Corn,Neib,Elm,P1,P2,P3)
Implicit None
!===========================================================================================
Intent(In)::Dim,Elm,P1,P2,P3
Intent(InOut)::Corn,Neib,NC

Integer::Dim,NC,P1,P2,P3,DT_Count,I,J,CommonCorners,Elm,E,N
Integer,Dimension(1:1000)::DT
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::isTriEdge,hasCorner
!===========================================================================================
if(hasCorner(Dim,Corn,Elm,P1) .And. hasCorner(Dim,Corn,Elm,P2) .And. hasCorner(Dim,Corn,Elm,P3)) then
	print *,'No Extra Triangle to Delete!!!'
else
!Part 1:
	DT_Count = 1
	DT(DT_Count) = Elm
	Call GetDeletingTris(Dim,Elm,Corn,Neib,DT,DT_Count,P1,P2,P3)

	if(DT_Count > 1) then
!Part 2:
		!-------------------------------- Defining new Tri cell --------------------------------
		NC = NC + 1
		Corn(NC,1) = P1 
		Corn(NC,2) = P2
		Corn(NC,3) = P3
		Corn(NC,4) = 0
		!------------------------- Initialization of neibours of new Tri -----------------------
		Neib(NC,1) = 0
		Neib(NC,2) = 0
		Neib(NC,3) = 0
		Neib(NC,4) = 0

		do I=1,DT_Count
			
			E = DT(I)
!Part 3:
			CommonCorners = 0
			do J=1,3
				if(Corn(E,J)==P1 .Or. Corn(E,J)==P2 .Or. Corn(E,J)==P3) then
					CommonCorners = CommonCorners + 1
				endif
			end do

			!----------------- Only Triangles with a common edge are considered ----------------
			if(CommonCorners == 2) then
!Part 4:
				if(isTriEdge(Corn(E,1),Corn(E,2),P1,P2,P3)) then
					
					N = Neib(E,3)

				elseif(isTriEdge(Corn(E,2),Corn(E,3),P1,P2,P3)) then

					N = Neib(E,1)

				elseif(isTriEdge(Corn(E,3),Corn(E,1),P1,P2,P3)) then
					
					N = Neib(E,2)

				endif

				do J=1,4
					if(Neib(N,J) == E) then
						Neib(N,J) = NC
						exit
					endif
				end do
				
				if(hasCorner(Dim,Corn,E,P1) .And. hasCorner(Dim,Corn,E,P2)) then

					Neib(NC,3) = N

				elseif(hasCorner(Dim,Corn,E,P2) .And. hasCorner(Dim,Corn,E,P3)) then
			
					Neib(NC,1) = N

				elseif(hasCorner(Dim,Corn,E,P3) .And. hasCorner(Dim,Corn,E,P1)) then
					
					Neib(NC,2) = N

				endif

			endif
!Part 5:
			Call DeleteCell(Dim,E,Corn)

		end do
		
	endif

endif
!===========================================================================================
End Subroutine TriangleFormation
!*********************************************************************************************
