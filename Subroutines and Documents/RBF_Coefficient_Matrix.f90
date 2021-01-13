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
!// Date: Feb., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Coefficient_Matrix(Dim,X,Y,IBP,NBP,A)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,IBP,NBP,X,Y
 Intent(InOut)::A

 Integer::Dim,I,J,P1,P2,NBP
 Real(8)::Dx,Dy,DL,Temp
 Integer,Dimension(1:Dim)::IBP
 Real(8),Dimension(1:Dim)::X,Y
 Real(8),Dimension(1:NBP+3,1:NBP+3)::A
!********************************************************************************************* 
!Part 1:
 Do I=1,NBP
    Do J=I+1,NBP

      !Part 2:
	   P1 = IBP(I)
	   P2 = IBP(J)

       Dx = X(P2)-X(P1)
       Dy = Y(P2)-Y(P1)
       DL = Dsqrt( Dx*Dx + Dy*Dy )

      !Part 3:
       Call RBF_Function(DL,Temp)
       A(I,J)=Temp 

    End Do
 End Do

!Part 4:
 Do I=1,NBP

   !Part 5:
    DL=0.0
    
   !Part 6:
	Call RBF_Function(DL,Temp)
    A(I,I)=Temp       

 End Do

!Part 7:
 Do I=1,NBP 
     A(I,NBP+1) = 1.    
     A(I,NBP+2) = X( IBP(I) )
     A(I,NBP+3) = Y( IBP(I) )
 End Do

!Part 8:
 A(NBP+1,NBP+1) = 0.0    
 A(NBP+1,NBP+2) = 0.0
 A(NBP+1,NBP+3) = 0.0
  
 A(NBP+2,NBP+2) = 0.0    
 A(NBP+2,NBP+3) = 0.0

 A(NBP+3,NBP+3) = 0.0

!Part 9:
 Do I=1,NBP+3
   Do J=I+1,NBP+3
      A(J,I)=A(I,J) 
   End Do
 End Do
!*********************************************************************************************
 End
!###########################################################################################