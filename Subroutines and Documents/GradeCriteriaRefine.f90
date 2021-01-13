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
 Subroutine GradeCriteriaRefine(Dim,NF,NBP,IBP,IDS,X,Y,CoRBF)
 Implicit None
!*********************************************************************************************
 Intent(In    ):: Dim,NF,NBP,IBP,IDS,X,Y
 Intent(Out )::CoRBF

 Integer::Dim,I,NF,NBP,P1,P2,NE,N
 Real(8)::X1,Y1,X2,Y2,EgeLen
 Integer,Dimension(Dim)::IBP
 Real(8),Dimension(Dim)::X,Y,Siz,CoRBF
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:NBP)::b
 Real(8),Dimension(1:NBP,1:NBP)::A
!*********************************************************************************************
!Part 1:
 Siz(:) = 0.0

!Part 2:
 Do I=1,NF
 
!Part 3:
 IF(IDS(2,I)/=0)Cycle

	  !Part 4:
       P1 = IDS(3,I)
       P2 = IDS(4,I)

       X1 = X(P1) ; Y1 = Y(P1)
       X2 = X(P2) ; Y2 = Y(P2)

	  !Part 5:
       EgeLen = sqrt(3.)*0.25*( (X2-X1)**2 + (Y2-Y1)**2 )

	  ! Part 6:
       Siz(P1) = Siz(P1) + EgeLen  !sqrt(3.0)/4.*EgeLen
       Siz(P2) = Siz(P2) + EgeLen  !sqrt(3.0)/4.*EgeLen

 End Do
 
!Part 7:
 Siz(:) = Siz(:)*0.5

!Part 8: 
 Do I=1,NBP
    N = IBP(I)
    b(I)=Siz(N)
 End Do

!Part 9:
 Call RBF_Coefficient_MatrixV2(Dim,NBP,IBP,X,Y,A)

!Part 10:
 Call solve_lu(NBP,A,b,CoRBF)

!*********************************************************************************************
    End
!###########################################################################################

