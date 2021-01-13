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
Recursive Subroutine GetDeletingTris(Dim,element,Corn,Neib,DT,DT_Count,P1,P2,P3)
Implicit None
!===========================================================================================
Intent(In)::Dim,element,Corn,Neib,P1,P2,P3
Intent(Inout)::DT,DT_Count

Integer::Dim,element,DT_Count,P1,P2,P3
Integer,Dimension(1:1000)::DT
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::IsInTheList,isTriEdge
!===========================================================================================
!Part 1:
if(.Not. isTriEdge(Corn(element,1),Corn(element,2),P1,P2,P3)) then

	!---------------------- Adding and checking other adjacent elements --------------------
	if(.Not. IsInTheList(DT,DT_Count,Neib(element,3))) then
		DT_Count = DT_Count + 1
		DT(DT_Count) = Neib(element,3)
		Call GetDeletingTris(Dim,Neib(element,3),Corn,Neib,DT,DT_Count,P1,P2,P3)
	endif

endif
!Part 2:
if(.Not. isTriEdge(Corn(element,2),Corn(element,3),P1,P2,P3)) then

	!---------------------- Adding and checking other adjacent elements --------------------
	if(.Not. IsInTheList(DT,DT_Count,Neib(element,1))) then
		DT_Count = DT_Count + 1
		DT(DT_Count) = Neib(element,1)
		Call GetDeletingTris(Dim,Neib(element,1),Corn,Neib,DT,DT_Count,P1,P2,P3)
	endif
endif
!Part 3:
if(.Not. isTriEdge(Corn(element,3),Corn(element,1),P1,P2,P3)) then

	!---------------------- Adding and checking other adjacent elements --------------------
	if(.Not. IsInTheList(DT,DT_Count,Neib(element,2))) then
		DT_Count = DT_Count + 1
		DT(DT_Count) = Neib(element,2)
		Call GetDeletingTris(Dim,Neib(element,2),Corn,Neib,DT,DT_Count,P1,P2,P3)
	endif

endif
!===========================================================================================
End Subroutine GetDeletingTris 
!*********************************************************************************************
