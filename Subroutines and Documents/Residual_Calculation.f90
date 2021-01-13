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
 Subroutine Residual_Calculation(Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,&
                                 GM,U0,V0,P0,R0,C0,WNP1,Residual)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,Neq,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2,IDS,NX,NY,DA,GM,U0,V0,P0,R0,C0,WNP1
 Intent(Out  )::Residual

 Integer::Dim,Neq,I,J,NC,NF,NF1,NF2,NFW1,NFW2,NFF1,NFF2
 Real(8)::GM,U0,V0,P0,R0,C0,U,V
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:Dim)::NX,NY,DA,P
 Real(8),Dimension(1:Neq,1:Dim)::WNP1,Con,Residual
 Real(8),Dimension(1:Neq+1,1:Dim)::WB
!*********************************************************************************************
!Part 1: 
 Do J=1,NC
	U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)
    P(J) = (GM-1)*(WNP1(4,J)-0.5*WNP1(1,J)*(U*U+V*V))
 End Do

!Part 2:
 Call BC_Wall(Dim,NFW1,NFW2,IDS,GM,P,WB)
 Call BC_Riemann(Dim,NFF1,NFF2,NX,NY,DA,IDS,GM,U0,V0,P0,R0,C0,WNP1,P,WB)
		
!Part 3:
 Call ConMeanFlow_AUSM(Dim,NC,NF1,NF2,NF,IDS,NX,NY,DA,GM,WNP1,WB,P,Con)
		
!Part 4:
 Do I=1,NC 
    Do J=1,4 
    Residual(J,I) = Con(J,I) 
    End Do    
 End Do
!*********************************************************************************************
 End
!###########################################################################################