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
!// Developed by: H. Morad Tabrizi, Mechanical Eng., Amirkabir University of Technology    //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Module LSQ																				 

Implicit None
!===============================
Public       :: STARTUP, INCLUD, REGCF, TOLSET, SING, SS, COV, PARTIAL_CORR, &
                VMOVE, REORDR, HDIAG, VARPRD, BKSUB2
!===============================
Private      :: INV
!===============================
Integer, Save, Public                     :: NOBS, NCOL, R_DIM
Integer, Allocatable, Save, Dimension(:), Public  :: VORDER, ROW_PTR
Logical, Save, Public                     :: INITIALIZED = .false., &
                                             TOL_SET = .false.,  &
                                             RSS_SET = .false.
!===============================
Integer, Parameter, Public      :: DP = Selected_Real_Kind(10,70)
Real (kind=DP), Allocatable, Save, Dimension(:), Public :: D, RHS, R, TOL, RSS
Real (kind=DP), Save, Private   :: ZERO = 0.0_DP, ONE = 1.0_DP, VSMALL
Real (kind=DP), Save, Public    :: SSERR, TOLY
!===============================
contains
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:........sUBRYTIN                           //!
!// Date         : May/1/2016                                                            //!
!// Developed by : H. Moradtabrizi, Iran, Tehran, OpenFlows@chmail.ir                    //!
!// Version      : V1                                                                    //!
!//                                                                                      //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.Marketxde.ir                       //!
!// It May Be xpied, Modified, And Redistributed For Non-xmmercial Use.                  //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************

Subroutine STARTUP(NVAR, FIT_CONST)

Integer, Intent(in)  :: NVAR
Logical, Intent(in)  :: FIT_CONST

Integer              :: I

VSMALL = 10.0_DP * tiny(ZERO)

NOBS = 0
If (FIT_CONST) Then
  NCOL = NVAR + 1
Else
  NCOL = NVAR
End If

If (INITIALIZED) Then
  Deallocate(D, RHS, R, TOL, RSS, VORDER, ROW_PTR)
End If
R_DIM = NCOL * (NCOL - 1)/2
Allocate( D(NCOL), RHS(NCOL), R(R_DIM), TOL(NCOL), RSS(NCOL), VORDER(NCOL),  &
          ROW_PTR(NCOL) )

D = ZERO
RHS = ZERO
R = ZERO
SSERR = ZERO

If (FIT_CONST) Then
  Do I = 1, NCOL
    VORDER(I) = I-1
  End Do
Else
  Do I = 1, NCOL
    VORDER(I) = I
  End Do
End If


ROW_PTR(1) = 1
Do I = 2, NCOL-1
  ROW_PTR(I) = ROW_PTR(I-1) + NCOL - I + 1
End Do
ROW_PTR(NCOL) = 0

INITIALIZED = .true.
TOL_SET = .false.
RSS_SET = .false.

Return
End Subroutine STARTUP

Subroutine INCLUD(WEIGHT, XROW, YELEM)

Real (kind=DP), Intent(in)                   :: WEIGHT, YELEM
Real (kind=DP), Dimension(:), Intent(in out) :: XROW


Integer         :: I, K, NEXTR
Real (kind=DP)  :: W, Y, XI, DI, WXI, DPI, CBAR, SBAR, XK

NOBS = NOBS + 1
W = WEIGHT
Y = YELEM
RSS_SET = .false.
NEXTR = 1
Do I = 1, NCOL

  If (ABS(W) < VSMALL) Then
    Return
  End If
  XI = XROW(I)
  if (ABS(XI) < VSMALL) Then
    NEXTR = NEXTR + NCOL - I
  Else
    DI = D(I)
    WXI = W * XI
    DPI = DI + WXI*XI
    CBAR = DI / DPI
    SBAR = WXI / DPI
    W = CBAR * W
    D(I) = DPI
    Do K = I+1, NCOL
      XK = XROW(K)
      XROW(K) = XK - XI * R(NEXTR)
      R(NEXTR) = CBAR * R(NEXTR) + SBAR * XK
      NEXTR = NEXTR + 1
    End Do
    XK = Y
    Y = XK - XI * RHS(I)
    RHS(I) = CBAR * RHS(I) + SBAR * XK
  End If
End Do 

SSERR = SSERR + W * Y * Y

Return
End Subroutine INCLUD

Subroutine REGCF(BETA, NREQ, IFAULT)


Integer, Intent(In)                        :: NREQ
Integer, Intent(Out)                       :: IFAULT
Real (kind=DP), Dimension(:), Intent(out)  :: BETA


Integer     :: I, J, NEXTR


IFAULT = 0
If (NREQ < 1 .or. NREQ > NCOL) Then
  IFAULT = IFAULT + 4
End If
If (IFAULT /= 0) Then
  Return
End If

If (.not. TOL_SET) Then
  Call TOLSET()
End If

Do I = NREQ, 1, -1
  If (SQRT(D(I)) < TOL(I)) Then
    BETA(I) = ZERO
    D(I) = ZERO
    IFAULT = -I
  Else
    BETA(I) = RHS(I)
    NEXTR = ROW_PTR(I)
    Do J = I+1, NREQ
      BETA(I) = BETA(I) - R(NEXTR) * BETA(J)
      NEXTR = NEXTR + 1
    End Do 
  End If
End Do 

Return

End Subroutine REGCF

Subroutine TOLSET(EPS)

Real (kind=DP), Intent(in), Optional :: EPS


Integer          :: COL, ROW, POS
Real (kind=DP)   :: EPS1, TOTAL
Real (kind=DP), Dimension(NCOL)  :: WORK
Real (kind=DP), Parameter        :: TEN = 10.0_DP


If (Present(EPS)) Then
  EPS1 = Max(ABS(EPS), TEN * Epsilon(TEN))
Else
  EPS1 = TEN * Epsilon(TEN)
End If

WORK = SQRT(D)
Do COL = 1, NCOL
  POS = COL - 1
  TOTAL = WORK(COL)
  Do ROW = 1, COL-1
    TOTAL = TOTAL + ABS(R(POS)) * WORK(ROW)
    POS = POS + NCOL - ROW - 1
  End Do
  TOL(COL) = EPS1 * TOTAL
End Do

TOL_SET = .true.
Return
End Subroutine TOLSET

Subroutine SING(LINDEP, IFAULT)

Integer, Intent(Out)                :: IFAULT
Logical, Dimension(:), Intent(Out)  :: LINDEP


Real (kind=DP)   :: TEMP, Y, WEIGHT
Integer          :: COL, POS, ROW, POS2
Real (kind=DP), Dimension(NCOL)  :: X, WORK

IFAULT = 0

WORK = SQRT(D)
If (.not. TOL_SET) Then
  Call TOLSET()
End If

Do COL = 1, NCOL
  TEMP = TOL(COL)
  POS = COL - 1
  Do ROW = 1, COL-1
    POS = POS + NCOL - ROW - 1
  End Do


  LINDEP(COL) = .false.
  If (WORK(COL) <= TEMP) Then
    LINDEP(COL) = .true.
    IFAULT = IFAULT - 1
    If (COL < NCOL) Then
      POS2 = POS + NCOL - COL + 1
      X = ZERO
      X(COL+1:NCOL) = R(POS+1:POS2-1)
      Y = RHS(COL)
      WEIGHT = D(COL)
      R(POS+1:POS2-1) = ZERO
      D(COL) = ZERO
      RHS(COL) = ZERO
      Call INCLUD(WEIGHT, X, Y)
                                                          
      NOBS = NOBS - 1
    Else
      SSERR = SSERR + D(COL) * RHS(COL)**2
    End If 
  End If 
End Do 
Return
End Subroutine SING

Subroutine SS()

Integer          :: I
Real (kind=DP)   :: TOTAL

TOTAL = SSERR
RSS(NCOL) = SSERR
Do I = NCOL, 2, -1
  TOTAL = TOTAL + D(I) * RHS(I)**2
  RSS(I-1) = TOTAL
End Do

RSS_SET = .true.
Return
End Subroutine SS

Subroutine COV(NREQ, VAR, COVMAT, DIMCOV, STERR, IFAULT)


Integer, Intent(In)                        :: NREQ, DIMCOV
Integer, Intent(Out)                       :: IFAULT
Real (kind=DP), intent(Out)                :: VAR
Real (kind=DP), dimension(:), intent(Out)  :: COVMAT, STERR


Integer           :: DIM_RINV, POS, ROW, START, POS2, COL, POS1, K
Real (kind=DP)    :: TOTAL
Real (kind=DP), Allocatable, Dimension(:)  :: RINV


If (DIMCOV < NREQ*(NREQ+1)/2) Then
  IFAULT = 1
  Return
End If


IFAULT = 0
Do ROW = 1, NREQ
  If (ABS(D(ROW)) < VSMALL) Then
    IFAULT = -ROW
  End If
End Do
If (IFAULT /= 0) Then
  Return
End If


If (NOBS > NREQ) Then
  If (.not. RSS_SET) Then
    Call SS()
  End If
  VAR = RSS(NREQ) / (NOBS - NREQ)
Else
  IFAULT = 2
  Return
End If

DIM_RINV = NREQ*(NREQ-1)/2
Allocate ( RINV(DIM_RINV) )

Call INV(NREQ, RINV)
POS = 1
START = 1
Do ROW = 1, NREQ
  POS2 = START
  Do COL = ROW, NREQ
    POS1 = START + COL - ROW
    If (ROW == COL) Then
      TOTAL = ONE / D(COL)
    Else
      TOTAL = RINV(POS1-1) / D(COL)
    End If
    Do K = COL+1, NREQ
      TOTAL = TOTAL + RINV(POS1) * RINV(POS2) / D(K)
      POS1 = POS1 + 1
      POS2 = POS2 + 1
    End Do
    COVMAT(POS) = TOTAL * VAR
    If (ROW == COL) Then
      STERR(ROW) = SQRT(COVMAT(POS))
    End If
    POS = POS + 1
  End Do
  START = START + NREQ - ROW
End Do

Deallocate(RINV)
Return
End Subroutine COV

Subroutine INV(NREQ, RINV)

Integer, Intent(In)                        :: NREQ
Real (kind=DP), Dimension(:), Intent(Out)  :: RINV

Integer          :: POS, ROW, COL, START, K, POS1, POS2
Real (kind=DP)   :: TOTAL


POS = NREQ * (NREQ-1)/2
Do ROW = NREQ-1, 1, -1
  START = ROW_PTR(ROW)
  Do COL = NREQ, ROW+1, -1
    POS1 = START
    POS2 = POS
    TOTAL = ZERO
    Do K = ROW+1, COL-1
      POS2 = POS2 + NREQ - K
      TOTAL = TOTAL - R(POS1) * RINV(POS2)
      POS1 = POS1 + 1
    End Do 
    RINV(POS) = TOTAL - R(POS1)
    POS = POS - 1
  End Do 
End Do 

Return
End Subroutine INV

Subroutine PARTIAL_CORR(INVAR, CORMAT, DIMC, YCORR, IFAULT)


Integer, Intent(In)                        :: INVAR, DIMC
Integer, Intent(Out)                       :: IFAULT
Real (kind=DP), Dimension(:), Intent(Out)  :: CORMAT, YCORR


Integer          :: BASE_POS, POS, ROW, COL, COL1, COL2, POS1, POS2
Real (kind=DP)   :: SUMXX, SUMXY, SUMYY
Real (kind=DP), Dimension(INVAR+1:NCOL)  :: RMS, WORK


IFAULT = 0
If (INVAR < 0 .or. INVAR > NCOL-1) Then
  IFAULT = IFAULT + 4
End If
If (DIMC < (NCOL-INVAR)*(NCOL-INVAR-1)/2) Then
  IFAULT = IFAULT + 8
End If
If (IFAULT /= 0) Then
  Return
End If


BASE_POS = INVAR*NCOL - (INVAR+1)*(INVAR+2)/2


If (D(INVAR+1) > ZERO) Then
  RMS(INVAR+1) = ONE / SQRT(D(INVAR+1))
End If
Do COL = INVAR+2, NCOL
  POS = BASE_POS + COL
  SUMXX = D(COL)
  Do ROW = INVAR+1, COL-1
    SUMXX = SUMXX + D(ROW) * R(POS)**2
    POS = POS + NCOL - ROW - 1
  End Do 
  If (SUMXX > ZERO) Then
    RMS(COL) = ONE / SQRT(SUMXX)
  Else
    RMS(COL) = ZERO
    IFAULT = -COL
  End If
End Do


SUMYY = SSERR
Do ROW = INVAR+1, NCOL
  SUMYY = SUMYY + D(ROW) * RHS(ROW)**2
End Do 
If (SUMYY > ZERO) Then
  SUMYY = ONE / SQRT(SUMYY)
End If

POS = 1
Do COL1 = INVAR+1, NCOL
  SUMXY = ZERO
  WORK(COL1+1:NCOL) = ZERO
  POS1 = BASE_POS + COL1
  Do ROW = INVAR+1, COL1-1
    POS2 = POS1 + 1
    Do COL2 = COL1+1, NCOL
      WORK(COL2) = WORK(COL2) + D(ROW) * R(POS1) * R(POS2)
      POS2 = POS2 + 1
    End Do
    SUMXY = SUMXY + D(ROW) * R(POS1) * RHS(ROW)
    POS1 = POS1 + NCOL - ROW - 1
  End Do


  POS2 = POS1 + 1
  Do COL2 = COL1+1, NCOL
    WORK(COL2) = WORK(COL2) + D(COL1) * R(POS2)
    POS2 = POS2 + 1
    CORMAT(POS) = WORK(COL2) * RMS(COL1) * RMS(COL2)
    POS = POS + 1
  End Do 
  SUMXY = SUMXY + D(COL1) * RHS(COL1)
  YCORR(COL1) = SUMXY * RMS(COL1) * SUMYY
End Do 

YCORR(1:INVAR) = ZERO

Return
End Subroutine PARTIAL_CORR

Subroutine VMOVE(FROM, DEST, IFAULT)


Integer, Intent(In)    :: FROM, DEST
Integer, Intent(Out)   :: IFAULT


Real (kind=DP)   :: D1, D2, X, D1NEW, D2NEW, CBAR, SBAR, Y
Integer          :: M, FIRST, LAST, INC, M1, M2, MP1, COL, POS, ROW

IFAULT = 0
If (FROM < 1 .or. FROM > NCOL) Then
  IFAULT = IFAULT + 4
End If
If (DEST < 1 .or. DEST > NCOL) Then
  IFAULT = IFAULT + 8
End If
If (IFAULT /= 0) Then
  Return
End If

If (FROM == DEST) Then
  Return
End If

If (.not. RSS_SET) Then
  Call SS()
End If

If (FROM < DEST) Then
  FIRST = FROM
  LAST = DEST - 1
  INC = 1
Else
  FIRST = FROM - 1
  LAST = DEST
  INC = -1
End If

Do M = FIRST, LAST, INC


  M1 = ROW_PTR(M)
  M2 = ROW_PTR(M+1)
  MP1 = M + 1
  D1 = D(M)
  D2 = D(MP1)


 If (D1 >= VSMALL .or. D2 >= VSMALL) Then
    X = R(M1)
    If (ABS(X) * SQRT(D1) < TOL(MP1)) Then
      X = ZERO
    End If
    If (D1 < VSMALL .or. ABS(X) < VSMALL) Then
      D(M) = D2
      D(MP1) = D1
      R(M1) = ZERO
      Do COL = M+2, NCOL
        M1 = M1 + 1
        X = R(M1)
        R(M1) = R(M2)
        R(M2) = X
        M2 = M2 + 1
      End Do
      X = RHS(M)
      RHS(M) = RHS(MP1)
      RHS(MP1) = X
    Else If (D2 < VSMALL) Then
      D(M) = D1 * X**2
      R(M1) = ONE / X
      R(M1+1:M1+NCOL-M-1) = R(M1+1:M1+NCOL-M-1) / X
      RHS(M) = RHS(M) / X
    Else


      D1NEW = D2 + D1*X**2
      CBAR = D2 / D1NEW
      SBAR = X * D1 / D1NEW
      D2NEW = D1 * CBAR
      D(M) = D1NEW
      D(MP1) = D2NEW
      R(M1) = SBAR
      Do COL = M+2, NCOL
        M1 = M1 + 1
        Y = R(M1)
        R(M1) = CBAR*R(M2) + SBAR*Y
        R(M2) = Y - X*R(M2)
        M2 = M2 + 1
      End Do
      Y = RHS(M)
      RHS(M) = CBAR*RHS(MP1) + SBAR*Y
      RHS(MP1) = Y - X*RHS(MP1)
    End If
  End If

  POS = M
  Do ROW = 1, M-1
    X = R(POS)
    R(POS) = R(POS-1)
    R(POS-1) = X
    POS = POS + NCOL - ROW - 1
  End Do 

  M1 = VORDER(M)
  VORDER(M) = VORDER(MP1)
  VORDER(MP1) = M1
  X = TOL(M)
  TOL(M) = TOL(MP1)
  TOL(MP1) = X
  RSS(M) = RSS(MP1) + D(MP1) * RHS(MP1)**2
End Do

Return
End Subroutine VMOVE

Subroutine REORDR(LIST, N, POS1, IFAULT)


Integer, Intent(in)                :: N, POS1
Integer, Dimension(:), Intent(in)  :: LIST
Integer, Intent(out)               :: IFAULT


Integer     :: NEXT, I, L, J
Logical     :: FOUND


IFAULT = 0
If (N < 1 .or. N > NCOL+1-POS1) Then
  IFAULT = IFAULT + 4
End If
If (IFAULT /= 0) Then
  Return
End If

NEXT = POS1
Do I = POS1, NCOL
  L = VORDER(I)
  FOUND = .false.
  Do J = 1, N
    If (L == LIST(J)) Then
      FOUND = .true.
      Exit
    End If
  End Do

  If (FOUND) Then
    If (I > NEXT) Then
      Call VMOVE(I, NEXT, IFAULT)
    End If
    NEXT = NEXT + 1
  End If
End Do

If (NEXT >= N+POS1) Then
  Return
End If

IFAULT = 8

Return
End Subroutine REORDR

Subroutine HDIAG(XROW, NREQ, HII, IFAULT)

Integer, Intent(In)                       :: NREQ
Integer, Intent(Out)                      :: IFAULT
Real (kind=DP), Dimension(:), Intent(In)  :: XROW
Real (kind=DP), Intent(Out)               :: HII


Integer         :: COL, ROW, POS
Real (kind=DP)  :: TOTAL
Real (kind=DP), Dimension(NCOL)  :: WK


IFAULT = 0
If (NREQ > NCOL) Then
  IFAULT = IFAULT + 4
End If
If (IFAULT /= 0) Then
  Return
End If


HII = ZERO
Do COL = 1, NREQ
  If (SQRT(D(COL)) <= TOL(COL)) Then
    WK(COL) = ZERO
  Else
    POS = COL - 1
    TOTAL = XROW(COL)
    Do ROW = 1, COL-1
      TOTAL = TOTAL - WK(ROW)*R(POS)
      POS = POS + NCOL - ROW - 1
    End Do 
    WK(COL) = TOTAL
    HII = HII + TOTAL**2 / D(COL)
  End If
End Do 

Return
End Subroutine HDIAG

Subroutine VARPRD(X, NREQ, FN_VAL)


Integer, Intent(In)                       :: NREQ
Real (kind=DP), Dimension(:), Intent(In)  :: X
Real (kind=DP), Intent(Out)               :: FN_VAL

Integer         :: IFAULT, ROW
Real (kind=DP)  :: VAR
Real (kind=DP), Dimension(NREQ) :: WK

FN_VAL = ZERO
IFAULT = 0
If (NREQ < 1 .or. NREQ > NCOL) Then
  IFAULT = IFAULT + 4
End If
If (NOBS <= NREQ) Then
  IFAULT = IFAULT + 8
End If
If (IFAULT /= 0) Then
  Write(unit=*, fmt="(1x, a, i4)") "Error in function VARPRD: ifault =", IFAULT
  Return
End If

VAR = SSERR / (NOBS - NREQ)

Call BKSUB2(X, WK, NREQ)
Do ROW = 1, NREQ
  If(D(ROW) > TOL(ROW)) Then
    FN_VAL = FN_VAL + WK(ROW)**2 / D(ROW)
  End If
End Do

FN_VAL = FN_VAL * VAR

Return
End Subroutine VARPRD

Subroutine BKSUB2(X, B, NREQ)

Integer, Intent(In)                        :: NREQ
Real (kind=DP), Dimension(:), Intent(In)   :: X
Real (kind=DP), Dimension(:), Intent(Out)  :: B

Integer           :: POS, ROW, COL
Real (kind=DP)   :: TEMP

Do ROW = 1, NREQ
  POS = ROW - 1
  TEMP = X(ROW)
  Do COL = 1, ROW-1
    TEMP = TEMP - R(POS)*B(COL)
    POS = POS + NCOL - COL - 1
  End Do
  B(ROW) = TEMP
End Do

Return
End Subroutine BKSUB2

End Module LSQ
!###########################################################################################
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!//////////////////////////////////////////////////////////////////////////////////////////!
!// Description:Polynomial fitting on a Function curve                                   //!                                                                                      //!
!//                                                                                      //!
!// Date: August,1,2016                                                                     //!
!// Developed by: H.Moradtabrizi, Iran, Tehran, h_mtabrizi@ut.ac.ir 					 //!
!// Version: V1                                                                          //!
!// Doc ID:                                                                              //!
!//                                                                                      //!
!// The Program Is Available Through The Website: www.MarketCode.ir                      //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                 //!
!//////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine Poly_Fit (N,X,Y,M,Pol_Coeff)
Use LSQ
Implicit None
!********************************************************************************************* 
Intent(In   )::N,X,Y,M
Intent(Out  )::Pol_Coeff
!===============================

Integer:: I,J,N,M,Ier
Real(8):: X(1:N),Y(1:N),Pol_Coeff(0:M),Var, Xrow(0:M+1),Wt = 1.0_dp
Logical:: Fit_const = .TRUE.      
!*********************************************************************************************

Call startup(M, Fit_const)
Do I = 2, N-1
  Xrow(0) = 1.0_dp
  Do J = 1, M
    Xrow(J) = X(I) * Xrow(J-1)
  End Do
  Call includ(Wt, Xrow, Y(I))
End Do

Call ss()
Var = rss(M+1) / (N - M - 1)



Call regcf(Pol_Coeff, M+1, Ier)


!###########################################################################################
 End 
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS