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
Subroutine getNextCorner(Dim,Corn,Elm,Node,P,Pi)
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
if(Corn(Elm,4) /= 0) then !---------------- Elm is a Quadrilateral ---------------------

    if(index == 4) then !------------ Finding Next corner to Node in Elm ---------------

	    Pi = 1            !--------- Index of Next corner after Node in Elm -----------
	    P = Corn(Elm,1)   !--------------- Next corner after Node in Elm --------------

    else

	    Pi = index+1      !---------- Index of Next corner after Node in Elm ----------
	    P = Corn(Elm,Pi)  !---------------- Next corner after Node in Elm -------------

    endif !---------------- End of Finding Next corner to Node in Elm ------------------

else !----------------------------------- Elm is a Triangle ----------------------------
!Part 3:
    if(index == 3) then !------------ Finding Next corner to Node in Elm ---------------

	    Pi = 1            !--------- Index of Next corner after Node in Elm -----------
	    P = Corn(Elm,1)   !--------------- Next corner after Node in Elm --------------

    else

	    Pi = index+1      !---------- Index of Next corner after Node in Elm ----------
	    P = Corn(Elm,Pi)  !---------------- Next corner after Node in Elm -------------

    endif !---------------- End of Finding Next corner to Node in Elm ------------------
    
endif

!===========================================================================================
End Subroutine getNextCorner
!*********************************************************************************************
