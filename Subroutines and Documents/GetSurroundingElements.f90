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
!// Date: April, 01, 2017                                                                  //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine GetSurroundingElements(Dim,Corn,Neib,P,NE,TElms,QElms,TEC,QEC)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,P,NE
Intent(Out)::TElms,QElms,TEC,QEC

Integer::Dim,TEC,QEC,P,NE,PreE,CurE,N1,N2,Next
Integer,Dimension(1:1000)::TElms,QElms
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::isBoundaryElement
!===========================================================================================
TEC = 0
QEC = 0
isBoundaryElement = .False.
!Part 1:
if(Corn(NE,4) == 0) then
	TEC = TEC + 1
	TElms(TEC) = NE
else
	QEC = QEC + 1
	QElms(QEC) = NE
endif

Call GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,NE,N1,N2)

if(N1 == N2) then !---------------------- Case When Only Two Element are Surrounding P ----------------------
!Part 2:
	if(N1 /= 0) then

		if(Corn(N1,4) == 0) then
			TEC = TEC + 1
			TElms(TEC) = N1
		else
			QEC = QEC + 1
			QElms(QEC) = N1
		endif

	endif

elseif(N1 == 0 .Or. N2 == 0) then
!Part 3:
	if(N1 /= 0) then

		PreE = NE
		CurE = N1

		do
			
			if(Corn(CurE,4) == 0) then
				TEC = TEC + 1
				TElms(TEC) = CurE
			else
				QEC = QEC + 1
				QElms(QEC) = CurE
			endif
			
			Call GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,CurE,N1,N2)
			
			if(N1 /= PreE) then
				if(N1 /= 0) then
					PreE = CurE
					CurE = N1
				else
					exit
				endif
			else
				if(N2 /= 0) then
					PreE = CurE
					CurE = N2
				else
					exit
				endif
			endif
			
		end do

	else !------------------------ N1 = 0 And N2 /= 0 ------------------------------

		PreE = NE
		CurE = N2

		do
			
			if(Corn(CurE,4) == 0) then
				TEC = TEC + 1
				TElms(TEC) = CurE
			else
				QEC = QEC + 1
				QElms(QEC) = CurE
			endif
			
			Call GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,CurE,N1,N2)
			
			if(N1 /= PreE) then
				if(N1 /= 0) then
					PreE = CurE
					CurE = N1
				else
					exit
				endif
			else
				if(N2 /= 0) then
					PreE = CurE
					CurE = N2
				else
					exit
				endif
			endif
			
		end do

	endif

else
!Part 4:
	PreE = NE
	CurE = N1 !----------------------- N1 or N2 Can be set here!!!!! -----------------------
	Next = N2

	do
		if(Corn(CurE,4) == 0) then
			TEC = TEC + 1
			TElms(TEC) = CurE
		else
			QEC = QEC + 1
			QElms(QEC) = CurE
		endif
		
		Call GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,CurE,N1,N2)
		
		if(N1 /= PreE) then
			if(N1 /= NE .And. N1 /= 0) then
				PreE = CurE
				CurE = N1
			else
				if(N1 == 0) then
					isBoundaryElement = .True.
				endif
				exit
			endif
		else
			if(N2 /= NE .And. N2 /= 0) then
				PreE = CurE
				CurE = N2
			else
				if(N2 == 0) then
					isBoundaryElement = .True.
				endif
				exit
			endif

		endif

	end do

	if(isBoundaryElement) then
		
		PreE = NE
		CurE = Next 

		do
			if(Corn(CurE,4) == 0) then
				TEC = TEC + 1
				TElms(TEC) = CurE
			else
				QEC = QEC + 1
				QElms(QEC) = CurE
			endif
			
			Call GetTwoNeibourSharingCorner(Dim,Corn,Neib,P,CurE,N1,N2)
			
			if(N1 /= PreE) then
				if(N1 /= 0) then
					PreE = CurE
					CurE = N1
				else
					exit
				endif
			else
				if(N2 /= 0) then
					PreE = CurE
					CurE = N2
				else
					exit
				endif

			endif

		end do

	endif

endif
!===========================================================================================
End Subroutine GetSurroundingElements 
!*********************************************************************************************
