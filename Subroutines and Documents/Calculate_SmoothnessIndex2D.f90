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
!// Developed by: M. Mohammadi, Mechanical Eng., Amirkabir University of Technology        //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Calculate_SmoothnessIndex2D(Dim,NC,X,Y,Corn,Neib,Smoothness)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,X,Y,Corn,Neib
 Intent(Out  )::Smoothness

 Integer::Dim,I,J,NC,P1,P2,P3,N1,N2,N3
 Real(8)::SumArea,X21,X31,Y21,Y31
 Integer,Dimension(1:4,1:Dim)::Corn,Neib
 Real(8),Dimension(1:Dim)::X,Y,Area,Smoothness
!*********************************************************************************************
!Part 1: 
 DO I=1,NC

    P1 = Corn(1,I)
    P2 = Corn(2,I)
    P3 = Corn(3,I)

    X21=X(P2)-X(P1)
    X31=X(P3)-X(P1)
    Y21=Y(P2)-Y(P1)
    Y31=Y(P3)-Y(P1)

    Area(I) = ABS(X21*Y31-X31*Y21) * 0.5
 ENDDO

!Part 1: 
 DO I=1,NC

   !Part 3:        
    DO J=1,3
		N1 = Neib(1,I)  
		N2 = Neib(2,I)  
		N3 = Neib(3,I)

       IF( N1==0 )THEN
        SumArea = ( Area(N2) + Area(N3) ) / 2 
	   ELSEIF( N2==0 )THEN
        SumArea = ( Area(N1) + Area(N3) ) / 2 
	   ELSEIF( N3==0 )THEN
        SumArea = ( Area(N1) + Area(N2) ) / 2 
	   ELSE
		SumArea = ( Area(N1) + Area(N2) + Area(N3) )/3
       ENDIF

    ENDDO

    Smoothness(I) = Area(I) / SumArea
 ENDDO
!*********************************************************************************************
 END