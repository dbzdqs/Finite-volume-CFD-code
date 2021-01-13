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
Function ElementInverted(Dim,Corn,X,Y,TElms,QElms,TEC,QEC,TEF,QEF)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,TElms,QElms,TEC,QEC

Integer::Dim,QEC,TEC,I
Integer,Dimension(1:1000)::TElms,QElms,TEF,QEF,TEF1,QEF1
Integer,Dimension(1:Dim,1:4)::Corn
Logical::ElementInverted
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:
ElementInverted = .False.
Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElms,QEC,TElms,TEC,TEF1,QEF1)

!Part 2:
do I=1,TEC
	if(TEF(I)/=TEF1(I)) then
		ElementInverted = .True.
		exit
	endif
end do

!Part 3:
if(.Not. ElementInverted) then
	do I=1,QEC
		if(QEF(I)/=QEF1(I)) then
			ElementInverted = .True.
			exit
		endif 
	end do
endif
!===========================================================================================
End Function ElementInverted
!*********************************************************************************************
