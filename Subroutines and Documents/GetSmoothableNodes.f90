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
!// Date: April, 01, 2017                                                                  //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine GetSmoothableNodes(Dim,NC,Corn,SmoothNode,PC)
Implicit None
!===========================================================================================
Intent(In)::Dim,NC,Corn
Intent(InOut)::SmoothNode,PC

Integer::Dim,NC,PC,I,J,Elm
Integer,Dimension(1:Dim,1:2)::SmoothNode
Integer,Dimension(1:Dim,1:4)::Corn
!===========================================================================================
!Part 1:
do I=1,NC

	if(Corn(I,1)/=Corn(I,2)) then

		Elm = I

		do J=1,4
			
            if(Corn(Elm,J) /= 0) then
                
			    Call AddToSmoothList(Dim,Corn(Elm,J),Elm,SmoothNode,PC)
            
            endif

		end do

	endif

end do
!===========================================================================================
End Subroutine GetSmoothableNodes
!*********************************************************************************************
