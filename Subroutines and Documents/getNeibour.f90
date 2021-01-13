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
Function getNeibour(Dim,Corn,Neib,P1,P2,E)
Implicit None
!===========================================================================================
Intent(In)::Dim,Corn,Neib,P1,P2,E

Integer::Dim,P1,P2,E,I,index1,index2,getNeibour
Integer,Dimension(1:Dim,1:4)::Corn,Neib
!===========================================================================================
if(E /= 0) then
!Part 1:
    if(Corn(E,4) /= 0) then
        
        do I=1,4

	        if(Corn(E,I) == P1) index1 = I
	        if(Corn(E,I) == P2) index2 = I

        end do

        if(index2 > index1) then

	        if(index2 - index1 == 1) then

		        getNeibour = Neib(E,index1)

	        else

		        getNeibour = Neib(E,index2)

	        endif

        else

	        if(index1 - index2 == 1) then

		        getNeibour = Neib(E,index2)

	        else

		        getNeibour = Neib(E,index1)

	        endif

        endif
    
    else
!Part 2:        
        do I=1,3
            
            if(Corn(E,I) /= P1 .And. Corn(E,I) /= P2) then
                index1 = I
                exit
            endif
            
        end do
        
        getNeibour = Neib(E,index1)
        
    endif

endif
!===========================================================================================
End Function getNeibour
!*********************************************************************************************
