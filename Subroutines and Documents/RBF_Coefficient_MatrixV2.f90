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
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Coefficient_MatrixV2(Dim,NBP,IBP,X,Y,A)
 Implicit None
!********************************************************************************************* 
 Intent(In   )::Dim,NBP,IBP,X,Y
 Intent(Out  )::A

 Integer::Dim,I,J,NBP,P1,P2
 Real(8)::Dx,Dy,DL,Temp
 Real(8),dimension(1:Dim)::X,Y
 Real(8),dimension(1:NBP,1:NBP)::A
 Integer,Dimension(Dim)::IBP 
!********************************************************************************************* 
!Part 1:
 Do I=1,NBP
    Do J=I+1,NBP

       P1 = IBP(I)
	   P2 = IBP(J)
       
      !Part 2:
       Dx = X(P2)-X(P1)
       Dy = Y(P2)-Y(P1)
       DL = Dsqrt( Dx*Dx + Dy*Dy )

      !Part 3:
       Call RBF_FunctionV2(DL,Temp)
       A(I,J)=Temp 

    End Do
 End Do

!Part 4:
 Do I=1,NBP

   !Part 5:
    DL=0.0
    
   !Part 6:
	Call RBF_FunctionV2(DL,Temp)
    A(I,I)=Temp       

 End Do

!Part 9:
 Do I=1,NBP
   Do J=I+1,NBP
      A(J,I)=A(I,J) 
   End Do
 End Do
 
!*********************************************************************************************
 End 
!###########################################################################################
