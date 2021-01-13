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
!// Date: Aug., 30, 2015                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function CurrentLevelCompleted(Dim,FrontEdges,Fronts,States,currLevel)
Implicit None
!===========================================================================================
Intent(In)::Dim,FrontEdges,Fronts,States,currLevel

Integer,Parameter::Processed = -1
Integer,Parameter::Level=3

Integer::Dim,Fronts,currLevel,I
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::FrontEdges
Logical::CurrentLevelCompleted
!===========================================================================================
!Part 1:
CurrentLevelCompleted = .True.

do I=1,Fronts
	if(States(I) /= Processed) then
		if(FrontEdges(I,Level) == currLevel) then
			CurrentLevelCompleted = .False.
			exit
		endif
	endif
end do
!===========================================================================================
End Function CurrentLevelCompleted 
!*********************************************************************************************
