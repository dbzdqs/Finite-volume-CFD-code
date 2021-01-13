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
Subroutine GetNeibouringPoints(Dim,Corn,Neib,QElms,QEC,NPList,NPC,V,Elm)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,QElms,QEC,V,Elm
Intent(InOut)::NPList,NPC

Integer::Dim,V,Elm,NPC,TEC,QEC,A,E,P,I,index,getNeibour
Integer,Dimension(1:1000)::QElms,NPList
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent
!===========================================================================================

NPC = 0 !------------------- Neibouring Points (of V) Counter ------------------------------

!Part 1:

Call getNextCorner(Dim,Corn,Elm,V,A,index)


Call addToList(NPList,NPC,A)

!Part 2:

P = A
E = Elm

do  
	
	do I=1,4
         
		if(areAdjacent(Dim,Corn,V,Corn(E,I),E) .And. Corn(E,I) /= P) then

			P = Corn(E,I)
			exit

        endif

	end do


	E = getNeibour(Dim,Corn,Neib,V,P,E)

	Call addToList(NPList,NPC,P)

	if(P == A) exit

end do

!===========================================================================================
End Subroutine GetNeibouringPoints
!*********************************************************************************************
