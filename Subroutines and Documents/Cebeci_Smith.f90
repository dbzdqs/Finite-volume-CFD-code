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
 Subroutine Cebeci_Smith(Dim,NC,NF1,NF2,NF,IDS,NX,NY,Xc,A,MR,NFW1,NFW2,DW,INW,Rinf,P,WB,&
                         WNP1,Mu,NNearstCell,INearstCell,Mut)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,Xc,A,MR,NFW1,NFW2,DW,INW,Rinf,P,WB,WNP1,Mu,&
                NNearstCell,INearstCell
 Intent(Out  )::Mut

 Integer::Dim,J,I,K,KK,P1,P2,NC,ME,NE,NF1,NF2,NF,II,NFW1,NFW2,JJ,I1,I2
 Real(8)::DX,UDX,Xr,RM,MR,Yn,Muw,Rw,Uw,Vw,Tw,Xm,XL,Vt,Rex,Delta,Yp,Ap,L,Li,L0,Mu0,ut,U,V,&
          DY,DL,Yk,TAUW,UNY,VNX,PNX,Rinf,AA,YMAX,UC_MAX,UC_MIN,UDIF,UC,RJ,E1,E2,&
		  FWAKE,BB,CC,FMAX,OMEGA,DEL,E,YN1,U_C,XJ,YJ,Xi,Yi,M_C,M,Dmin,Dis,PP,MUT_MAX
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::NX,NY,Xc,A,Mu,Mut,DW,DUDY,DVDX,DPDX,MUTO,MUTI,FKLEB,UE,RE,DELSTAR,P,YM
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
 Integer,Dimension(NFW1:NFW2,1:12000)::INearstCell
 Integer,Dimension(NFW1:NFW2)::NNearstCell
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    DUDY(I) = 0.0
    DVDX(I) = 0.0
    DPDX(I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U  = WB(2,I)/WB(1,I)
	V  = WB(3,I)/WB(1,I)
    PP = WB(5,I)

   !Part 5:
    DUDY(ME) = DUDY(ME) + U *NY(I)
    DVDX(ME) = DVDX(ME) + V *NX(I)
    DPDX(ME) = DPDX(ME) + PP*NX(I)
 End Do

!Part 6:
 DO I=NF1+1,NF2

   !Part 7:
    ME = IDS(1,I)
    NE = IDS(2,I)

   !Part 8:
    U   = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )
    V   = 0.5*( WNP1(3,ME)/WNP1(1,ME) + WNP1(3,NE)/WNP1(1,NE) )
    PP  = 0.5*( P(ME) + P(NE) )
    
    UNY = U  * NY(I)
    VNX = V  * NX(I)
    PNX = PP * NX(I)
    
    DUDY(ME) = DUDY(ME) + UNY
    DUDY(NE) = DUDY(NE) - UNY

    DVDX(ME) = DVDX(ME) + VNX
    DVDX(NE) = DVDX(NE) - VNX

    DPDX(ME) = DPDX(ME) + PNX
    DPDX(NE) = DPDX(NE) - PNX

 End Do

!Part 9:
 DO I=1,NC
    DUDY(I) = DUDY(I)/A(I)
    DVDX(I) = DVDX(I)/A(I)
    DPDX(I) = DPDX(I)/A(I)
 End Do

!Part 10:
 Xr=0.0

!Part 11:calculate the velocity on the boundary layer edge = UE  
 DO K=NFW1+1,NFW2

    !Part 12:
     E1=10000.0
	 

    !Part 13:
     ME = IDS(1,K)

    !Part 14:
     XL = Dabs( Xc(ME)-Xr )

    !Part 15:
     Rex = XL * Rinf

    !Part 16:
     Delta = DABS ( 0.16*XL*Rex**(-1./7.) )

   !Part 17:
    DO J=1,NC

      !Part 18:
       I=INW(J)
       Yn = DW(J)

      !Part 19:
       IF(I/=K) CYCLE
    
      !Part 20:
       IF(DABS(Delta-Yn) < E1  ) THEN
        U     = WNP1(2,J)/WNP1(1,J)
        V     = WNP1(3,J)/WNP1(1,J)
      
        E1 = DABS(Delta-Yn)
 
        UE(K)= DSQRT(U*U + V*V)
        RE(K)  = WNP1(1,J)
    
       END IF 
	  
    END DO

 END DO

!Part 21:calculate the displacement thickness = DELSTAR   
 DO K=NFW1+1,NFW2

   !Part 22: 
    DELSTAR(K) = 0.0

   !Part 23:
    ME = IDS(1,K)

   !Part 24:
    XL = Dabs( Xc(ME)-Xr )

   !Part 25:
  	Rex = XL * Rinf

   !Part 26:
    Delta = DABS ( 0.16*XL*Rex**(-1./7.) )
 
   !Part 27: 
    DO J=1,NC

      !Part 28:
       I=INW(J)
       Yn=DW(J)
     
      !Part 29:
       IF (I/=K .and. Yn>Delta) CYCLE

      !Part 30:
       RJ  = WNP1(1,J)
       U   = WNP1(2,J)/WNP1(1,J)
       V   = WNP1(3,J)/WNP1(1,J)

      !Part 31:
       U_C = DSQRT(U*U + V*V)
	   
      !Part 32:    
       DEL = DABS ( (1- ( (RJ*U_C) /(UE(K)*RE(K) ) )) * Yn)
       DELSTAR(K) = DELSTAR(K) + DEL  

    END DO

   !Part 33:   
    IF(delstar(k)==0.0 .OR. delstar(k)>DELTA ) delstar(k)=0.5*DELTA
 END DO    

!Part 34:
 RM = 1/MR

!Part 35: 
 Do J=1,NC

   !Part 36:
    Yn = DW(J)
    I  = INW(J)

   !Part 37:
    ME = IDS(1,I)

   !Part 38:
    Rw  = WNP1(1,ME) 
    Muw = Mu(ME)
    TAUW= Muw*( DUDY(ME) )
	
   !Part 39:
    IF(TAUW==0.0) THEN
	  Ap = 26.0
	 ELSEIF(DPDX(J)==0) THEN
	  Ap = 26.0
	 ELSE
      Ap = 26.0/DSQRT(1+DABS(Yn*DPDX(J)/TAUW))
    END IF
    
   !Part 40:
    XL = Dabs( Xc(ME)-Xr )

   !Part 41:
    Rex = XL * Rinf

   !Part 42:
    Delta = 0.16*XL*Rex**(-1./7.) 

   !Part 43:
    ut = Dsqrt( Dabs(TAUW)/Rw )
    Yp = Rw*ut*Yn/Muw  * Dsqrt(RM)

   !Part 44:
    L = 0.4 * Yn * ( 1-exp(-Yp/Ap) )

   !Part 45:
    OMEGA = Dsqrt( DUDY(J)**2 + DVDX(J)**2)
    MutI(J) = RM *  L*L * WNP1(1,J) * OMEGA

         
   !Part 46:
    FKLEB(J) = 1/(1+5.5*((Yn/Delta)**6))
    
   !Part 47:
    MutO(J) = 0.0168 * RM * FKLEB(J) * WNP1(1,J) *  UE(I) * DELSTAR(I) 
  
 End Do

!Part 48:
 YM(:)=0.0

!Part 49:
 DO K=NFW1+1,NFW2

    DO JJ=2,NNearstCell(K)
       I1=INearstCell(K,JJ-1)
       I2=INearstCell(K,JJ)
        
	   IF( ( (MutO(i1) - MutI(i1)) * (MutO(i2) - MutI(i2)) )<0. )then
	     YM(K) = DW(i1)
		 exit
	endif

	END DO

 END DO

!Part 50:
 DO J=1,NC

   !Part 51:
    Yn = DW(J)
    I  = INW(J)

   !Part 52:
    IF(Yn<YM(I)) THEN
	   Mut(J) = MutI(J)
    ELSE
	   Mut(J) = MutO(J)
    END IF


 END DO

!Part 53:
 DO II=NFW1+1,NFW2
    
   !Part 54:
    MUT_MAX=0.0

  !Part 55:
   DO J=1,NC

    I=INW(J)
    IF(I/=II) cycle
    IF(MUT_MAX < Mut(J)) MUT_MAX = Mut(J)

   END DO

  !Part 56:
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

