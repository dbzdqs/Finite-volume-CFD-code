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
 Subroutine Baldwin_Lomax(Dim,NC,NF1,NF2,NF,NFW1,NFW2,IDS,NX,NY,A,INW,DW,NNearstCell,INearstCell,MR,Mut0,WB,WNP1,Mu,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,NFW1,NFW2,IDS,NX,NY,DW,NNearstCell,INearstCell,MR,INW,Mut0,WB,WNP1,Mu
 Intent(Out  )::Mut

 Integer::Dim,J,I,K,II,P1,P2,NC,ME,NE,NF1,NF2,NF,NFW1,NFW2,I1,JJ,I2,I3,I4
 Real(8)::DX,UDX,RM,MR,Yn,Muw,Rw,Uw,Vw,Tw,Ap,L,Mu0,ut,U,V,TAUW,UNY,VNX,UC_MAX,UC_MIN,UC,&
		  FWAKE,M,M_C,Fmaxi,Ymaxi,UDIFi,FKLEB,Fy_MAX,OMEGA,MUT_MAX,Xj,Yj,Xi,Yi,DY,E,&
		  Mut0,Xr,XL,Rex,DEL,Delta
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::NX,NY,A,Mu,Mut,MutI,MutO,Fmax,Ymax,Yp,DUDY,DVDX,DW,UDIF,DUY,YM,UE,Fy
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(NFW1:NFW2,1:12000)::INearstCell
 Integer,Dimension(NFW1:NFW2)::NNearstCell
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
 RM = 1/MR
 Ap = 26.0

!Part 11:
 Do II=NFW1+1,NFW2

   !Part 12:
    UC_MIN = 10000.0
    UC_MAX = 0.0
    Fy_MAX = 0.0

   !Part 13:
    ME = IDS(1,II)
    Rw  = WNP1(1,ME) 
    Muw = Mu(ME)
    TAUW= Muw * DUDY(ME) 
    ut = Dsqrt( Dabs(TAUW)/Rw )

   !Part 14:
    Do J=1,NC

      !Part 15:
       Yn = DW(J)
       I  = INW(J)

       IF(I/=II)Cycle

        Yp(J) = Rw*ut*Yn/Muw  * Dsqrt(RM)
 
       !Part 16:
        OMEGA =  Dabs( DUDY(J)   - DVDX(J) )
        Fy(J) = Yn * OMEGA * ( 1-exp(-Yp(J)/Ap) )

      !Part 17:
       U = WNP1(2,J)/WNP1(1,J)
       V = WNP1(3,J)/WNP1(1,J)

       UC = Dsqrt( U*U + V*V )
 
       IF(UC > UC_MAX) UC_MAX = UC
       IF(UC < UC_MIN) UC_MIN = UC

      !Part 18:
       IF(Fy(J) > Fy_MAX)Then
	    Fy_MAX = Fy(J)
        YMAXi = Yn
       Endif

    End Do

   !Part 19:
    UDIF(II) = Dabs( UC_MAX  - UC_MIN ) 
    Ymax(II)=Ymaxi
    Fmax(II)=Fy_MAX
 
 End Do

!Part 20:
 DO K=NFW1+1,NFW2
	
    DO JJ=3,NNearstCell(K)
      I1=INearstCell(K,JJ-2)
      I2=INearstCell(K,JJ-1)
      I3=INearstCell(K,JJ)
	   
	   IF( ( (Fy(i2) - Fy(i1)) * (Fy(i3) - Fy(i2)) )<0. )then
	    Fmax(K) = Fy(I2)
		Ymax(K) = DW(I2)
                
       !Part 21:
		I4=INearstCell(K,JJ+4)
		IF(Fy(I4)< Fy(I2)) exit
       Endif
 		
     END DO

 END DO

!Part 22:
 Do J=1,NC

   !Part 23:
    Yn = DW(J)
    I  = INW(J)

   !Part 24:
    Fmaxi=Fmax(I)
    Ymaxi=Ymax(I)
    UDIFi=UDIF(I)

   !Part 25:
    L = 0.4 * Yn * ( 1-exp(-Yp(J)/Ap) )

   !Part 26:
    OMEGA = Dabs( DUDY(J)  - DVDX(J)  )

   !Part 27:
    MutI(J) = RM *  L*L * WNP1(1,J) * OMEGA

   !Part 28:
    FKLEB = 1/(1+5.5*((Yn*0.3/Ymaxi)**6))

   !Part 29:
    IF(FMAXi==0.0) THEN
     MutO(J)=0.0 !tozih dar file marbote
    ELSE
	 FWAKE = min( Fmaxi*Ymaxi , 0.25 * YMAXi * UDIFi * UDIFi / FMAXi )   
     MutO(J)= RM *   0.0168 * 1.6 * WNP1(1,J) * FWAKE * FKLEB
    END IF

 END DO

!Part 30:
 YM(:)=10000.0
 DO K=NFW1+1,NFW2

    DO JJ=2,NNearstCell(K)
       I1=INearstCell(K,JJ-1)
       I2=INearstCell(K,JJ)

	   IF( ( (MutO(i1) - MutI(i1)) * (MutO(i2) - MutI(i2)) )<0. )Then
	    YM(K) = DW(i1)
	    Exit
	   Endif

	END DO

 END DO

!Part 31:
 DO J=1,NC

   !Part 32:
    Yn = DW(J)
    I  = INW(J)

   !Part 33:
    IF(YN<YM(I)) THEN
     Mut(J) = MutI(J)
    ELSE
     Mut(J) = MutO(J)
    END IF

 END DO

!Part 34:
 DO II=NFW1+1,NFW2
    
   !Part 35:
    MUT_MAX=0.0

   !Part 36:
    DO J=1,NC
       I=INW(J)
       IF(I/=II) cycle
       IF(MUT_MAX < Mut(J)) MUT_MAX = Mut(J)
    END DO

   !Part 37:
    IF(MUT_MAX < 8) THEN
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
