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
!// Developed by: A. Rezaii, Maritime eng., Amirkabir University of Technology             //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine RBF_Moving_Mesh(Dim,NBP,NP,IBP,X,Y,DelX,DelY)
 Implicit None
!********************************************************************************************* 
 Intent(In    ):: Dim,NBP,NP,IBP,X,Y
 Intent(InOut )::DelX,DelY

 Integer::Dim,I,N,NBP,NP  ,j
 Real(8),Dimension(Dim)::X,Y,DelX,DelY
 Integer,Dimension(1:Dim)::IBP
 Real(8),Dimension(1:NBP+3)::Cox,Coy,b_x,b_y
 Real(8),Dimension(1:NBP+3,1:NBP+3)::A
!********************************************************************************************* 
!Part 1:
 Do I=1,NBP
    N = IBP(I)

    b_x(I)=DelX(N)
    b_y(I)=DelY(N)
 End Do

 b_x(NBP+1)=0.0 
 b_x(NBP+2)=0.0 
 b_x(NBP+3)=0.0 

 b_y(NBP+1)=0.0 
 b_y(NBP+2)=0.0 
 b_y(NBP+3)=0.0 


!Part 2:
 Call RBF_Coefficient_Matrix(Dim,X,Y,IBP,NBP,A)

!Part 3:
!Call CG_Matrix_Solver(Dim,NBP,A,b_x,Cox)
!Call CG_Matrix_Solver(Dim,NBP,A,b_y,Coy)

!Solving by LU decomposition
 Call solve_lu(NBP+3,A,b_x,Cox)
 Call solve_lu(NBP+3,A,b_y,Coy)

!Part 4:
 Call RBF_Point_Dispalcement(Dim,IBP,NP,NBP,Cox,Coy,X,Y,DelX,DelY)
!*********************************************************************************************
 End
!###########################################################################################