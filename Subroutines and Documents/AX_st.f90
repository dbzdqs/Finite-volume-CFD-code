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
!// Developed by: S. Sheikhi, petrolium, Amirkabir university of Technology                //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine AX_st(Dim,n,nz_num,ia,ja,a,x,w)
Implicit None
!**********************************************************************************************
Intent(In)::Dim,N,NZ_NUM,IA,JA,A,X
Intent(Out  )::W

Integer::N,NZ_NUM,I,J,K,Dim
Real(8),Dimension(1:Dim)::A
Real(8),Dimension(1:4,1:Dim)::X,W
Integer,Dimension(1:Dim)::IA,JA
!**********************************************************************************************
!Part1:
 W(1:4,1:N) = 0.0D+00

Do K = 1, nz_num
  
  !Part2:
   I = IA(K)
   J = JA(K)

  !Part3:
   W(1,I) =W(1,I) + A(K) * X(1,J)
   W(2,I) =W(2,I) + A(K) * X(2,J)
   W(3,I) =W(3,I) + A(K) * X(3,J)
   W(4,I) =W(4,I) + A(K) * X(4,J)

 end Do
!*********************************************************************************************   
 End
!###########################################################################################