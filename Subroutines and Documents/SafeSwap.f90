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
Subroutine SafeSwap(Dim,Corn,Neib,X,Y,FrontEdges,Fronts,States,tri_down,tri_up)
Implicit None
!===========================================================================================
Intent(In)::Dim,Fronts,States,X,Y,tri_down,tri_up
Intent(Inout)::FrontEdges,Corn,Neib

Integer,Parameter::Non_Front=-1
Integer,Parameter::LeftVertex=1
Integer,Parameter::RightVertex=2
Integer,Parameter::Element=4

Integer::Dim,Fronts,tri_down,tri_up
Integer,Dimension(1:Dim)::States
Integer,Dimension(1:Dim,1:4)::Corn,Neib,FrontEdges
Logical::SwapIsSafe
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
if(SwapIsSafe(Dim,Corn,X,Y,tri_down,tri_up)) then

	print *,'tri_d1:',tri_down,Corn(tri_down,1),Corn(tri_down,2),Corn(tri_down,3)
	print *,'tri_u1:',tri_up,Corn(tri_up,1),Corn(tri_up,2),Corn(tri_up,3)
!Part 1:
	Call Swap(Dim,tri_down,tri_up,Corn,Neib)

	print *,'tri_d2:',tri_down,Corn(tri_down,1),Corn(tri_down,2),Corn(tri_down,3)
	print *,'tri_u2:',tri_up,Corn(tri_up,1),Corn(tri_up,2),Corn(tri_up,3)
!Part 2:
	Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,tri_down)
	Call UpdateFrontElements(Dim,Corn,FrontEdges,Fronts,States,tri_up)

endif
!===========================================================================================
End Subroutine SafeSwap
!*********************************************************************************************
