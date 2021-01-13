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
!// Developed by: *//*-+/01                        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Function areAdjacent(Dim,Corn,P1,P2,E)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,P1,P2,E

Integer::Dim,P1,P2,E,index1,index2,I
Integer,Dimension(1:Dim,1:4)::Corn
Logical::areAdjacent
!===========================================================================================
areAdjacent = .False.
!Part 1:
if(Corn(E,4) /= 0) then
    
    do I=1,4
	
	    if(Corn(E,I) == P1) then
		    index1 = I
	    endif

	    if(Corn(E,I) == P2) then
		    index2 = I
	    endif

    end do

    if(index1 == 1 .Or. index1 == 3) then
	    if(index2 == 4 .Or. index2 == 2) then
		    areAdjacent = .True.
	    endif
    elseif(index1 == 2 .Or. index1 == 4) then
	    if(index2 == 1 .Or. index2 == 3) then
		    areAdjacent = .True.
	    endif
    endif

else
!Part 2:    
    if(P1 /= P2) then
        
        areAdjacent = .True.
    
    else
    
        areAdjacent = .False.
        
    endif
    
endif
!===========================================================================================
End Function areAdjacent
!*********************************************************************************************
