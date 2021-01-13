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
Subroutine getPrevCorner(Dim,Corn,Elm,Node,P,Pi)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Elm,Node
Intent(Out)::P,Pi

Integer::Dim,Elm,Node,P,Pi,index,I
Integer,Dimension(1:Dim,1:4)::Corn
!===========================================================================================
!Part 1:
do I=1,4
	
	if(Corn(Elm,I) == Node) then

		index = I
		exit

	endif
	
end do
!Part 2:
if(index == 1) then !---------- Finding Previous corner to Node in Elm --------------

	Pi = 4 		      !-------- Index of Previous corner before Node in Elm --------
	P = Corn(Elm,4)   !------------ Previous corner before Node in Elm -------------

else

	Pi = index-1      !------- Index of Previous corner before Node in Elm ---------
	P = Corn(Elm,Pi)  !------------ Previous corner before Node in Elm -------------

endif !----------- End of Finding Previous corner to Node in Elm --------------------

!===========================================================================================
End Subroutine getPrevCorner
!*********************************************************************************************
