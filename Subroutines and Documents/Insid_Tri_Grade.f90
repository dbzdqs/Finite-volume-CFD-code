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
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Insid_Tri_Grade(Dim,NC,X,Y,SF,Corn,Grad,Cand_Tti)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,X,Y,SF,Corn
 Intent(Out  )::Cand_Tti
 Intent(InOut)::Grad

 Integer::Dim,J,P1,P2,P3,NC,Cand_Tti
 Real(8)::Max,X12,Y12,X13,Y13,X23,Y23,SF_C,Edge_12,Edge_13,Edge_23
 Integer,Dimension(1:Dim,1:4)::Corn
 Integer,Dimension(1:Dim)::Grad
 Real(8),Dimension(1:Dim)::X,Y,SF,Area
!*********************************************************************************************
!Part 1: 
 Cand_Tti=0 
 Max=0.0

!Part 2:
 Do J=1,NC
    If( Grad(J)==0 )Then 

    !Part 3:
     P1 = Corn(J,1)
     P2 = Corn(J,2)
	 P3 = Corn(J,3)

    !Part 4:
     X12 = X(P2) - X(P1)
	 Y12 = Y(P2) - Y(P1)

     X13 = X(P3) - X(P1)
	 Y13 = Y(P3) - Y(P1)

     X23 = X(P3) - X(P2)
	 Y23 = Y(P3) - Y(P2)

    !Part 5:
     Area(J) = Dabs( X12*Y13 - Y12*X13 ) / 2

    !Part 6:
     SF_C = ( SF(P1) + SF(P2) + SF(P3) ) / 3

    !Part 7:
     Edge_12=Dsqrt( X12*X12 + Y12*Y12 )
	 Edge_13=Dsqrt( X13*X13 + Y13*Y13 )
	 Edge_23=Dsqrt( X23*X23 + Y23*Y23 )

    !Part 8:
	 If( Area(J)<(SF_C*SF_C/2) .And. Edge_12<2*SF_C .And. Edge_13<2*SF_C .And. Edge_23<2*SF_C )Then
	  Grad(J) = 1
	 Else
	  Grad(J) =-1
	 End If

    Endif
 End Do

!Part 9: 
 Do J=1,NC

    If( Grad(J)==-1 .And. Area(J)>Max )Then 
     Max = Area(J)
	 Cand_Tti = J
	End If

 End Do

!*********************************************************************************************
 End 
!###########################################################################################
