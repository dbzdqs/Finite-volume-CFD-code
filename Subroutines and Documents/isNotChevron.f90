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
Function isNotChevron(Dim,Corn,X,Y,Q)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,X,Y,Q

Integer::Dim,A,B,C,D,Q
Integer,Dimension(1:Dim,1:4)::Corn
Logical::IntersectionOccur,isNotChevron
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
if(Q /= 0)then
!Part 1:    
    A = Corn(Q,1)
    B = Corn(Q,2)
    C = Corn(Q,3)
    D = Corn(Q,4)

    if(IntersectionOccur(Dim,A,C,B,D,X,Y)) then
	    isNotChevron = .True.
    else
	    isNotChevron = .False.
    endif

endif
!===========================================================================================
End Function isNotChevron
!*********************************************************************************************
