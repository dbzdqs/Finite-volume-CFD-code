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
!// Date: Feb., 20, 2016                                                                   //!
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
 Subroutine JacobianVector_Product(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,&
                                   GM,U0,V0,P0,R0,C0,DT,&
								   WNP1,z,Az) 
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,A,GM,U0,V0,P0,R0,C0,DT,WNP1,z
 Intent(Out  )::Az

 Integer::Dim,Neq,I,J,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2
 Real(8)::GM,U0,V0,P0,R0,C0,eps,norm_w,norm_z,Temp,Machin_eps
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,A,DA,DT
 Real(8),Dimension(1:Neq,1:Dim)::Az,z,WNP1,WNP1_Plus,WNP1_Minus,Rn_Plus,Rn_Minus
!********************************************************************************************* 
!Part 1:
 norm_z = 0.0
 norm_W = 0.0
  
 Do I=1,NC
    Do J=1,Neq
       norm_z = norm_z + z(J,I)*z(J,I)
       norm_W = norm_W + WNP1(J,I)*WNP1(J,I)
    End Do 
 End Do
 norm_z=Sqrt(norm_z)
 norm_W=Sqrt(norm_W)

!Part 2:
 Machin_eps = epsilon(0.0)
 IF( norm_z>Machin_eps ) eps=  Sqrt( (1.0+norm_W)*Machin_eps ) / norm_z
 
!Part 3:
 IF( norm_z<Machin_eps  )Then
  
  Az=0.0

 Else
    
 !Part 4:
  Do I=1,NC
     Do J=1,Neq
        WNP1_Plus(J,I)  = WNP1(J,I) + eps*z(J,I)
        WNP1_Minus(J,I) = WNP1(J,I) - eps*z(J,I)
     End Do 
  End Do
 
 !Part 5: 
  Call Residual_Calculation(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,GM,U0,V0,P0,R0,C0,WNP1_Plus ,Rn_Plus )
  Call Residual_Calculation(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,GM,U0,V0,P0,R0,C0,WNP1_Minus,Rn_Minus)
 
 !Part 6:
  Do I=1,NC

     Temp = 3.0*A(I)/DT(I)

     Do J=1,Neq
        Az(J,I) = Temp * z(J,I)  +  ( Rn_Plus(J,I) - Rn_Minus(J,I) ) / eps 
     End Do 

  End Do  
 
 End IF
!********************************************************************************************* 
 End
!###########################################################################################
