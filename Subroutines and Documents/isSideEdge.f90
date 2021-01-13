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
Function isSideEdge(newQuad,Ei)
Implicit None
!===========================================================================================
Intent(In)::newQuad,Ei

Integer::Dim,Fronts,I,P,Q
Integer,Dimension(1:2)::Ei
Integer,Dimension(1:4)::newQuad
Logical::isSideEdge
!===========================================================================================
isSideEdge = .False.
do I=1,4
!Part 1:
    if(I /= 4) then
        P = newQuad(I)
        Q = newQuad(I+1)
    else
        P = newQuad(I)
        Q = newQuad(1)    
    endif
!Part 2:
    if(Ei(1) == P .And. Ei(2) == Q) then
        isSideEdge = .True.
        exit
    elseif(Ei(1) == Q .And. Ei(2) == P) then
        isSideEdge = .True.
        exit
    endif
end do
!===========================================================================================
End Function isSideEdge
!*********************************************************************************************
