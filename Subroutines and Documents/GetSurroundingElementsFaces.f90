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
Subroutine GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms,QEC,TElms,TEC,TEF,QEF)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,QElms,QEC,TElms,TEC
Intent(Out)::TEF,QEF

Integer::Dim,QEC,TEC,E,I,J,C,TriSide
Integer,Dimension(1:4)::QSide
Integer,Dimension(1:1000)::QElms,TElms,TEF,QEF
Integer,Dimension(1:Dim,1:4)::Corn
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
do I=1,TEC
	E = TElms(I)
	TEF(I) = TriSide(Dim,Corn,X,Y,E) 
end do
!Part 2:
do I=1,QEC

	C = 0
	E = QElms(I)
	Call QuadSide(Dim,Corn,X,Y,E,QSide)
		
	do J=1,4
		if(QSide(J)==-1) C=C+1
	end do

	if(C>1) then
		QEF(I) = -1	
	else
		QEF(I) = 1
	endif

end do
!===========================================================================================
End Subroutine GetSurroundingElementsFaces
!*********************************************************************************************
