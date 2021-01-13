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
Subroutine CalcGamma(Dim,Corn,Neib,X,Y,V,Elm,Gdir_x,Gdir_y,GAMMA)
!===========================================================================================
Intent(In)::Dim,Corn,Neib,X,Y,V,Elm,Gdir_x,Gdir_y
Intent(Out)::GAMMA

Integer::Dim,V,Elm,C,Ci,P,E,I,getNeibour
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::areAdjacent
Real(8)::GAMMA,Gdir_x,Gdir_y,SE,GetNorm,temp,G_Norm
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================
!Part 1:

Call getNextCorner(Dim,Corn,Elm,V,C,Ci)
!-------------------- Finding shortest edge connected to 'V' corner ------------------------
SE = GetNorm(Dim,V,C,X,Y)
P = C
E = Elm

!Part 2:

do
            
	do I=1,4

		if(areAdjacent(Dim,Corn,V,Corn(E,I),E) .And. Corn(E,I) /= P) then
							
			P = Corn(E,I)
            E = getNeibour(Dim,Corn,Neib,V,P,E)
            temp = GetNorm(Dim,V,P,X,Y) 
            
            if(temp < SE) then
                SE = temp    
            endif
            
			exit

		endif

	end do

	if(P == C) exit

end do

!Part 3:

!----------------------------- Calculate norm of Gradient ----------------------------------

G_Norm = DSQRT(Gdir_x*Gdir_x + Gdir_y*Gdir_y)

GAMMA = (SE*0.001)/G_Norm

!===========================================================================================
End Subroutine CalcGamma 
!*********************************************************************************************
