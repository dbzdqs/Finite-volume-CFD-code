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
Subroutine setNeibour(Dim,Corn,Neib,ME,P1,P2,NE)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,ME,P1,P2,NE
Intent(Out)::Neib

Integer::Dim,ME,P1,P2,NE,I,index1,index2
Integer,Dimension(1:Dim,1:4)::Corn,Neib
!===========================================================================================
if(ME /= 0) then
!Part 1:    
    if(Corn(ME,4) /= 0) then

	    do I=1,4

		    if(Corn(ME,I) == P1) index1 = I
		    if(Corn(ME,I) == P2) index2 = I

	    end do

	    if(index1 > index2) then

		    if(index1 - index2 == 1) then

			    Neib(ME,index2) = NE

		    else

			    Neib(ME,index1) = NE

		    endif

	    else

		    if(index2 - index1 == 1) then

			    Neib(ME,index1) = NE

		    else

			    Neib(ME,index2) = NE

		    endif

	    endif

    else
!Part 2:
	    do I=1,3

		    if(Corn(ME,I) /= P1 .And. Corn(ME,I)/= P2) then
		
			    index1 = I
			    exit

		    endif

	    end do

	    Neib(ME,index1) = NE

    endif

endif
!===========================================================================================
End Subroutine setNeibour
!*********************************************************************************************
