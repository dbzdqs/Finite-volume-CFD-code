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
!// Date: Dec., 05, 2016                                                                   //!
!// Developed by: K. Moradian, Computer Science, Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine CollapseOperation(Dim,Corn,Neib,X,Y,NBE,BFP,D,A,B,ME,done)
Implicit None
!===========================================================================================
Intent(In)::Dim,NBE,BFP,ME
Intent(InOut)::Corn,Neib,X,Y,done

Integer::Dim,NBE,ME,I,J,A,Ai,B,Bi,C,Ci,D,Di,Nb,Nd,TEC,Nab,Nbc,Ncd,Nda,getNeibour
Integer,Dimension(1:1000)::TElms,QElmsB,QElmsD,TEF,QEF_B,QEF_D
Integer,Dimension(1:Dim,1:2)::BFP
Integer,Dimension(1:Dim,1:4)::Corn,Neib
Logical::done,ElementInverted,IsOnExteriorBoundary
Real(8)::x_B,y_B,x_D,y_D,x_value,y_value
Real(8),Dimension(1:Dim)::X,Y
!===========================================================================================

done = .False.

!Part 1:

if(.Not. IsOnExteriorBoundary(Dim,NBE,BFP,B) .And. .Not. IsOnExteriorBoundary(Dim,NBE,BFP,D)) then
    
    do I=1,4
	
	    if(Corn(ME,I) == A) Ai = I
	    if(Corn(ME,I) == B) Bi = I
	    if(Corn(ME,I) == D) Di = I

    end do

!Part 2:
    
    Call GetSurroundingElements(Dim,Corn,Neib,B,ME,TElms,QElmsB,TEC,Nb)
    Call GetSurroundingElements(Dim,Corn,Neib,D,ME,TElms,QElmsD,TEC,Nd)

    Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElmsB,Nb,TElms,TEC,TEF,QEF_B)
    Call GetSurroundingElementsFaces(Dim,Corn,X,Y,QElmsD,Nd,TElms,TEC,TEF,QEF_D)

!Part 3:
    
    x_B = X(B)
    y_B = Y(B)

    x_D = X(D)
    y_D = Y(D)
	
    x_value = (X(B) + X(D))/2
    y_value = (Y(B) + Y(D))/2

    X(B) = x_value
    Y(B) = y_value
    X(D) = x_value
    Y(D) = y_value

!Part 4:
    
    if(ElementInverted(Dim,Corn,X,Y,TElms,QElmsD,TEC,Nd,TEF,QEF_D) .Or. ElementInverted(Dim,Corn,X,Y,TElms,QElmsB,TEC,Nb,TEF,QEF_B)) then
	
	    X(B) = x_B
	    Y(B) = y_B
							
	    X(D) = x_D
	    Y(D) = y_D

    else

	    done = .True.

	    Call getOppoCorner(Dim,Corn,ME,A,C,Ci)

	    Nab = getNeibour(Dim,Corn,Neib,A,B,ME)
	    Nda = getNeibour(Dim,Corn,Neib,A,D,ME)
	
	    Nbc = getNeibour(Dim,Corn,Neib,B,C,ME)
	    Ncd = getNeibour(Dim,Corn,Neib,C,D,ME)

	    Call setNeibour(Dim,Corn,Neib,Nab,A,B,Nda)
	    Call setNeibour(Dim,Corn,Neib,Nda,D,A,Nab)
	    Call setNeibour(Dim,Corn,Neib,Nbc,B,C,Ncd)
	    Call setNeibour(Dim,Corn,Neib,Ncd,C,D,Nbc)

	    do I=1,Nd
		    do J=1,4
			    if(Corn(QElmsD(I),J) == D) then
				    Corn(QElmsD(I),J) = B	
			    endif
		    end do
	    end do

	    !-------------------- Deleting ME --------------------

	    do I=2,4
		    Corn(ME,I) = Corn(ME,1)
	    end do

	    !---------->>>>>>> Smoothing neibours >>>>>>>>>>>>>>>>

	    do I=1,4

		    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nab,I),Nab)

	    end do

	    do I=1,4

		    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nbc,I),Nbc)

	    end do

	    do I=1,4

		    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Ncd,I),Ncd)

	    end do

	    do I=1,4

		    Call ConstrainedLaplacianSmooth(Dim,NBE,BFP,Corn,Neib,X,Y,Corn(Nda,I),Nda)

	    end do

    endif

endif
!===========================================================================================
End Subroutine CollapseOperation 
!*********************************************************************************************
