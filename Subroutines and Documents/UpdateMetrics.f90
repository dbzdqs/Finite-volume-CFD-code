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
Subroutine UpdateMetrics(Dim,Corn,X,Y,QElms,Mu_min,E_min,min_index,QEC,PreDistortionMetrics,COINCIDENT_TOLERANCE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,QElms,QEC,COINCIDENT_TOLERANCE
Intent(InOut)::PreDistortionMetrics,Mu_min,E_min,min_index

Integer::Dim,I,min_index,E_min,QEC
Integer,Dimension(1:1000)::QElms
Integer,Dimension(1:Dim,1:4)::Corn
Real(8)::Mu_min,QuadDistortionMetric,COINCIDENT_TOLERANCE
Real(8),Dimension(1:100)::PreDistortionMetrics
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
do I=1,QEC

	PreDistortionMetrics(I) = QuadDistortionMetric(Dim,Corn,X,Y,QElms(I),COINCIDENT_TOLERANCE) 

end do

Mu_min = PreDistortionMetrics(1)
E_min = QElms(1)
min_index = 1

do I=2,QEC

	if(Mu_min > PreDistortionMetrics(I)) then
	
		Mu_min = PreDistortionMetrics(I)
		E_min = QElms(I)
		min_index = I

	endif

end do 
!===========================================================================================
End Subroutine UpdateMetrics
!*********************************************************************************************
