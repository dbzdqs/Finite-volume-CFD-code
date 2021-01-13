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
!// Date: Mar., 10, 2015                                                                   //!
!// Developed by: N. msnkre, Aerospace Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine TimSTP_Inviscid3D(Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,DA,CFLx,GM,P,WNP1,WB,DT)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,NZ,Vol,DA,CFLx,GM,P,WNP1,WB
 Intent(Out  )::DT

 Integer::Dim,I,NC,ME,NE,NF1,NF2,NF,j
 Real(8)::U,V,W,T1,R,R1,R2,C,DX,CFLx,GM
 Integer,Dimension(1:6,1:Dim)::IDS
 Real(8),Dimension(1:6,1:Dim)::WB
 Real(8),Dimension(1:5,1:Dim)::WNP1
 Real(8),Dimension(1:Dim)::DT,Vol,P,NX,NY,NZ,DA
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
    W = WB(4,I)/R

   !Part 5:
    C = SQRT( ABS( GM*WB(6,I)/R ) )

   !Part 6:
    DT(ME) = DT(ME) + ABS(U*NX(I) + V*NY(I) + W*NZ(I)) + C*DA(I)

 End Do

!Part 7:
 DO I=NF1+1,NF2
 
   !Part 8:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 9:
    R1 = WNP1(1,ME)
	R2 = WNP1(1,NE)
    U = 0.5*( WNP1(2,ME)/R1 + WNP1(2,NE)/R2 )
    V = 0.5*( WNP1(3,ME)/R1 + WNP1(3,NE)/R2 )
    W = 0.5*( WNP1(4,ME)/R1 + WNP1(4,NE)/R2 )

   !Part 10:
    C = SQRT( ABS( GM*(P(ME)+P(NE))/(R1+R2) ) )

   !Part 11:
    T1 = ABS(U*NX(I) + V*NY(I) + W*NZ(I)) + C*DA(I)
 
   !Part 12:
    DT(ME) = DT(ME) + T1
    DT(NE) = DT(NE) + T1
 End Do
    
!Part 13:
 DO I=1,NC
    DT(I) = CFLx*Vol(I)/DT(I)
 End Do
!*********************************************************************************************
 End
!###########################################################################################
