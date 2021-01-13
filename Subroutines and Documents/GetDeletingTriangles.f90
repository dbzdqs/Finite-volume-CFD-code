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
Recursive Subroutine GetDeletingTriangles(Dim,element,Corn,Neib,DT,DT_Count,newQuad,possible)
Implicit None
!===========================================================================================
Intent(In)::Dim,element,Corn,Neib,newQuad
Intent(Inout)::DT,DT_Count,possible

Integer::Dim,element,DT_Count,I
Integer,Dimension(1:4)::newQuad
Integer,Dimension(1:1000)::DT
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::isInserted,isQuadEdge,possible
!===========================================================================================
if(element /= 0) then
    
    if(Corn(element,4) == 0) then
!Part 1:    
        if(.not. (isQuadEdge(Corn(element,1),Corn(element,2),newQuad))) then

	        !---------------- Checking if target element has already been considered or not --------
	        isInserted = .false.
	        do I=1,DT_Count 
		        if(DT(I) == Neib(element,3)) then
			        isInserted = .true.
			        exit
		        endif
	        end do
	        !---------------------- Adding and checking other adjacent elements --------------------
	        if(.Not. isInserted) then
		        DT_Count = DT_Count + 1
		        DT(DT_Count) = Neib(element,3)
		        Call GetDeletingTriangles(Dim,Neib(element,3),Corn,Neib,DT,DT_Count,newQuad,possible)
	        endif
        endif
!Part 2:
        if(.not. (isQuadEdge(Corn(element,2),Corn(element,3),newQuad))) then

	        !---------------- Checking if target element has already been considered or not --------
	        isInserted = .false.
	        do I=1,DT_Count 
		        if(DT(I) == Neib(element,1)) then
			        isInserted = .true.
			        exit
		        endif
	        end do
	        !---------------------- Adding and checking other adjacent elements --------------------
	        if(.Not. isInserted) then
		        DT_Count = DT_Count + 1
		        DT(DT_Count) = Neib(element,1)
		        Call GetDeletingTriangles(Dim,Neib(element,1),Corn,Neib,DT,DT_Count,newQuad,possible)
	        endif
        endif
!Part 3:
        if(.not. (isQuadEdge(Corn(element,3),Corn(element,1),newQuad))) then

	        !---------------- Checking if target element has already been considered or not --------
	        isInserted = .false.
	        do I=1,DT_Count 
		        if(DT(I) == Neib(element,2)) then
			        isInserted = .true.
			        exit
		        endif
	        end do
	        !---------------------- Adding and checking other adjacent elements --------------------
	        if(.Not. isInserted) then
		        DT_Count = DT_Count + 1
		        DT(DT_Count) = Neib(element,2)
		        Call GetDeletingTriangles(Dim,Neib(element,2),Corn,Neib,DT,DT_Count,newQuad,possible)
	        endif

        endif
    
    else
        possible = .False.    
    endif

else
    possible = .False.    
endif
!===========================================================================================
End Subroutine GetDeletingTriangles
!*********************************************************************************************
