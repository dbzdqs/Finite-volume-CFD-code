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
!// Developed by: M. Keley, Mechanical Eng., Amirkabir University of Technology            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Thomas(Dim,NC,NF,NF1,NF2,NFW1,NFW2,IDS,NX,NY,A,INW,Mut0,Rinf,MR,Mu,WB,WNP1,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF,NF1,NF2,NFW1,NFW2,IDS,Rinf,NX,NY,A,WB,Mut0,INW,WNP1,MR,Mu
 Intent(Out  )::Mut

 Integer::Dim,I,K,ME,NE,NC,NF,NF1,NF2,J,II,NFW1,NFW2,P1,P2
 Real(8)::UNY,VNX,UC_MIN,UC_MAX,OME,UC,OMEGAC,RM,MR,Li,L0,L1,Delta,XL,Rex,Rinf,Yn,Xr,Ap,&
          Rw,ut,TAUW,Muw,Yp,U,V,Dmin,Dis,OMEGA,DY,DX,MUT_MAX,Mut0
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::NX,NY,A,Mut,DUDY,DVDX,Mu,L
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
!*********************************************************************************************
!Part 1:
 DO I=1,NC
    DUDY(I) = 0.0
    DVDX(I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)
    V = WB(3,I)/WB(1,I)

   !Part 5:
    DUDY(ME) = DUDY(ME) + U*NY(I)
    DVDX(ME) = DVDX(ME) + V*NX(I)

 End Do

!Part 6:
 DO I=NF1+1,NF2

   !Part 7:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 8:
    U   = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V   = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )

    UNY = U*NY(I)
	VNX = V*NX(I)

    DUDY(ME) = DUDY(ME) + UNY
    DUDY(NE) = DUDY(NE) - UNY

    DVDX(ME) = DVDX(ME) + VNX
    DVDX(NE) = DVDX(NE) - VNX

 End Do

!Part 9:
 DO I=1,NC
    DUDY(I) = DUDY(I)/A(I)
    DVDX(I) = DVDX(I)/A(I)
 End Do

!Part 10:
 Do II=NFW1+1,NFW2

!Part 11:
 UC_MIN = 10000.0
 UC_MAX = 0.0
 OME = 0.0
 OMEGAC= 0.0

!Part 12:
 Do J=1,NC

   !Part 13:
    I  = INW(J)

    IF(I/=II)Cycle
 
   !Part 14:
    OMEGA =  Dabs( DUDY(J)   - DVDX(J) )

   !Part 15:
    U = WNP1(2,J)/WNP1(1,J)
    V = WNP1(3,J)/WNP1(1,J)

    UC = Dsqrt( U*U + V*V )
 
   !Part 16: 
    IF(UC > UC_MAX) UC_MAX = UC
    IF(UC < UC_MIN) UC_MIN = UC

   !Part 17: 
    IF(OMEGAC < OMEGA ) OMEGAC=OMEGA

 End Do

!Part 18:
 L(II) = 0.13*(UC_MAX - UC_MIN ) / OMEGAC

 End Do

!part 19:  
 RM = 1/MR 
                                 
!Part 20:
 DO J=1,NC

    I=INW(J)

    OMEGA =  Dabs( DUDY(J)   - DVDX(J) )
    Mut(J) =  WNP1(1,J)* L(I) * L(I) * OMEGA * RM

 END DO

!Part 21:
 DO II=NFW1+1,NFW2

  !Part 22:
   MUT_MAX=0.0

 !Part 23:
  DO J=1,NC
   I=INW(J)
   IF(I/=II) cycle
   IF(MUT_MAX < Mut(J)) MUT_MAX = Mut(J)
  END DO

 !Part 24:
  IF(MUT_MAX < 8.0) THEN
      DO J=1,NC
       I=INW(J)
       IF(I/=II) cycle
       Mut(J) = 0.0
      END DO
  END IF

 END DO
!*********************************************************************************************
 End
!###########################################################################################

