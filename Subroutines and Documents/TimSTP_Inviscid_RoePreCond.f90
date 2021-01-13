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
 Subroutine TimSTP_Inviscid_RoePreCond(Dim,NC,NF1,NF2,NF,IDS,NX,NY,A,DA,CFLx,GM,P,WNP1,WB,Minf,DT)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,CFLx,IDS,A,P,GM,WNP1,Wb,Minf
 Intent(Out  )::DT

 Integer::Dim,I,NC,ME,NE,NF1,NF2,NF
 Real(8)::U,V,T,T1,R,R1,R2,C,C_h,DX,CFLx,GM,M,M2,Minf,U_h,kk
 Integer,Dimension(1:4,1:Dim)::IDS
 Real(8),Dimension(1:5,1:Dim)::Wb
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::DT,A,P,NX,NY,DA
!*********************************************************************************************	
!Part 1:
 Do I=1,NC
    DT(I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    R = WB(1,I)
    U = WB(2,I)/R
    V = WB(3,I)/R

   !Part 5:
    C = SQRT( ABS( GM*WB(5,I)/R ) )
	M=(U*U+V*V)/(C*C)

   !Part 6:
    M2=max(Minf*Minf,M)
    T=min(M2,1.0)

   !Part 7: 
	U_h=0.5*(1+T)*(U*NX(I)+V*NY(I))
	C_h=0.5*SQRT(4*C*C*T+(1-T)*(1-T)*(U*NX(I)+V*NY(I))*(U*NX(I)+V*NY(I)))

   !Part 8:
	T1 = ABS(U_h) + C_h*DA(I)
    DT(ME) = DT(ME) + T1

 End Do

!Part 9:
 DO I=NF1+1,NF2
 
   !Part 10:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 11:
    R1 = WNP1(1,ME)
	R2 = WNP1(1,NE)
    U = 0.5*( WNP1(2,ME)/R1 + WNP1(2,NE)/R2 )
    V = 0.5*( WNP1(3,ME)/R1 + WNP1(3,NE)/R2 )

   !Part 12:
    C = SQRT( ABS( GM*(P(ME)+P(NE))/(R1+R2) ) )
	M=(U*U+V*V)/(C*C)

   !Part 13:       
	M2=max(Minf*Minf,M)
    T=min(M2,1.0)

   !Part 14:
	U_h=0.5*(1+T)*(U*NX(I)+V*NY(I))
	C_h=0.5*SQRT(4*C*C*T+(1-T)*(1-T)*(U*NX(I)+V*NY(I))*(U*NX(I)+V*NY(I)))

   !Part 15:
	T1 = ABS(U_h) + C_h*DA(I)
	

   !Part 16:
    DT(ME) = DT(ME) + T1
    DT(NE) = DT(NE) + T1
 End Do

!Part 17:
 DO I=1,NC
    DT(I) = CFLx*A(I)/DT(I)
 End Do
!*********************************************************************************************
 End
!###########################################################################################