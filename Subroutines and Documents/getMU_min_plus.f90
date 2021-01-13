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
Subroutine getMU_min_plus(PreDistortionMetrics,QEC,Gradient,Gdir_x,Gdir_y)
Implicit None
!===========================================================================================
Intent(In)::PreDistortionMetrics,QEC,Gradient
Intent(Out)::Gdir_x,Gdir_y

Integer,Parameter::Gx = 1
Integer,Parameter::Gy = 2

Integer::QEC,I,J,index
Real(8)::Gdir_x,Gdir_y,key,G_norm,temp
Real(8),Dimension(1:100)::PreDistortionMetrics,Values
Real(8),Dimension(1:100,1:2)::Gradient
!===========================================================================================

!Part 1:

do I=1,QEC !--------------- Since PreDistortionMetrics Should'nt be Modified --------------

	Values(I) = PreDistortionMetrics(I) 	

end do

!---------- Sorting Distortion Metrics using 'INSERSTION SORT' (ascending order) ----------

do I=2,QEC

	key = Values(I)
	J = I-1

	do

        if(J < 1) then
            exit    
        elseif(Values(J)<key) then
            exit    
        endif
        
		Values(J+1) = Values(J)
		J = J-1

	end do

	Values(J+1) = key

end do

!Part 2:

do I=2,QEC

	temp = Values(I)

	do J=1,QEC

		if(temp == PreDistortionMetrics(J)) then

			index = J
			exit

		endif

	end do

	Gdir_x = Gradient(index,Gx)
	Gdir_y = Gradient(index,Gy)

	G_norm = DSQRT(Gdir_x*Gdir_x + Gdir_y*Gdir_y)
	
	if(G_norm > 0.1) exit	

end do

!===========================================================================================
End Subroutine getMU_min_plus
!*********************************************************************************************
