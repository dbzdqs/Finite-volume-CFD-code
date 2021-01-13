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
!// Supervisor:                                                                            //!
!// Chief Developer: N. msnkre, Aerospace eng. Amirkabir University of Technology          //!
!// Date: Feb., 10, 2018                                                                   //!
!// Supervisor:                                                                            //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
 Subroutine Prandtl_Escudier(Dim,NC,NF1,NF2,NF,IDS,NX,NY,Xc,A,INW,DW,MR,Rinf,Mut0,WB,WNP1,Mu,Mut,DUY,R0,U0,Mu0)
 Implicit None
!*********************************************************************************************
 Intent(In   )::Dim,NC,NF1,NF2,NF,IDS,NX,NY,Xc,A,INW,DW,MR,Rinf,Mut0,WB,WNP1,Mu,DUY,R0,U0,Mu0
 Intent(Out  )::Mut

 Integer::Dim,J,I,NearestW,P1,P2,NC,ME,NE,NF1,NF2,NF,ii
 Real(8)::U,DX,UDX,Xr,RM,MR,Yn,Muw,Rw,Uw,Vw,Tw,Xm,XL,Vt,Rex,Delta,Yp,Ap,L,Li,L0,Utau,&
          X21,X10,Y21,Y10,DY,DL,Ym,Yk,TAUW,UNY,Rinf,Mut0,R0,U0,Mu0
 Integer,Dimension(1:4,1:Dim)::IDS
 Integer,Dimension(1:Dim)::INW
 Real(8),Dimension(1:Dim)::NX,NY,Xc,A,Mu,Mut,DW,DUDY,DUY
 Real(8),Dimension(1:4,1:Dim)::WNP1
 Real(8),Dimension(1:5,1:Dim)::WB
!*********************************************************************************************	
!Part 1:
 DO I=1,NC
    DUDY(I) = 0.0
 End Do

!Part 2:
 DO I=NF2+1,NF

   !Part 3:
    ME = IDS(1,I)

   !Part 4:
    U = WB(2,I)/WB(1,I)

   !Part 5:
    DUDY(ME) = DUDY(ME) + U*NY(I)
    
 End Do

!Part 6:
 DO I=NF1+1,NF2

   !Part 7:
    ME = IDS(1,I)
	NE = IDS(2,I)

   !Part 8:
    U   = 0.5*( WNP1(2,ME)/WNP1(1,ME) + WNP1(2,NE)/WNP1(1,NE) )

    UNY = U*NY(I)

    DUDY(ME) = DUDY(ME) + UNY
    DUDY(NE) = DUDY(NE) - UNY

 End Do

!Part 9:
 DO I=1,NC
    DUDY(I) = DUDY(I)/A(I)
 End Do

!Part 10:
 Xr = 0.0
 RM = 1/MR
 Ap = 26.0

!Part 11:
 Do J=1,NC

   !Part 12:
	Yn       = DW(J)
    NearestW = INW(J)

   !Part 13:
    ME = IDS(1,NearestW)

   !Part 14:
	Rw  = WNP1(1,ME)
	Uw  = WNP1(2,ME)/Rw
	Vw  = WNP1(3,ME)/Rw 
    Muw = Mu(ME)
	TauW= Muw*   DUDY(ME) ! DUY(NearestW)  !

   !Part 15:
    XL = Dabs(Xc(ME)-Xr) 

   !Part 16:
	Rex = XL * Rinf

   !Part 17:
    !Delta = 0.16*XL*(Rex * RM)**(-1.0/7.0)
   Delta = 0.16*XL*(Rex)**(-1.0/7.0)

   Delta = 0.16* XL*( R0*U0*XL/Mu0 *RM) **(-1./7)


   !Part 18:
	L0 = 0.089 * Delta 

   !Part 19:
    Utau = Dsqrt( Dabs( TauW)/Rw )
    Yp   = Rw*Utau*Yn/Muw * Dsqrt(RM)

   !Part 20:
	Li = 0.42 * Yn * ( 1.0 - exp(-Yp/Ap) )

   !Part 21:
    L = Dmin1(Li,L0)

   !Part 22:
	Mut(J) = RM * L*L * WNP1(1,J) * Dabs( DUDY(J) )

   !Part 23:
	IF(Mut(J)<Mut0 ) Mut(J) = Mut0
 !If(NearestW==32840)print*,yn , L0
 End Do


    UNY=0.0

   !Part 6:
    Do I=1,NC

	   If(Mut(I)>UNY)then
	    UNY = Mut(J)
	    j   = I
	   Endif

	End do

 !print*,UNY,j
 !pause
!*********************************************************************************************
 End
!###########################################################################################
	  
    
