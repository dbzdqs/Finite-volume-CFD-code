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
!// Developed by: K. Safari, Mathmatical, Amirkabir university of Technology               //!
!//                                                                                        //!
!// The Program is Available Through the Website: www.DANS.ir                              //!
!// It May be Copied, Modified and Redistributed for Non-Commercial Use.                   //!
!//----------------------------------------------------------------------------------------//!
!// Duty:                                                                                  //!
!//                                                                                        //!
!////////////////////////////////////////////////////////////////////////////////////////////!
!*********************************************************************************************
Subroutine SVD(a,MDim,w,matu,u,matv,v,ierr)
 Implicit None
!*********************************************************************************************
 Integer::MDim   
 Real(8), INTENT(IN),DIMENSION(MDim,MDim) :: a
 Real(8),INTENT(OUT),DIMENSION(MDim):: w
 LOGICAL,INTENT(IN):: matu
 Real(8),INTENT(OUT),DIMENSION(MDim,MDim):: u
 LOGICAL,INTENT(IN):: matv
 Real(8),INTENT(OUT),DIMENSION(MDim,MDim):: v
 INTEGER,INTENT(OUT):: ierr
 Real(8) :: c,f,g,h
 INTEGER :: i,j,k,l, i1,k1
 INTEGER:: its
 INTEGER:: m,n
 Real(8),DIMENSION(SIZE(a,2)):: rv1
 Real(8):: s,x,y,z,scale,anorm
 Real(8),PARAMETER:: ZERO = 0.0, ONE=1.0, TWO=2.0
!*********************************************************************************************
!Part 1:
 w(:)   = 0
 u(:,:) = 0
 v(:,:) = 0
  
 m = SIZE(a,1)
 n = SIZE(a,2)
 IF (SIZE(w) < n) THEN
  ierr = -1
  RETURN
 END IF
  
 IF (matu .AND. (SIZE(u,1) < m .OR. SIZE(u,2) < n) ) THEN
  ierr = -2
  RETURN
 END IF
  
 IF (matv .AND. (SIZE(v,1) < n .OR. SIZE(v,2) < n) ) THEN
  ierr = -3
  RETURN
 END IF
    
 ierr = 0
 u(1:m,1:n) = a
 g     = ZERO
 scale = ZERO
 anorm = ZERO

!Part 2:
 DO  i = 1, n
    rv1(i) = scale * g
    g      = ZERO
    s      = ZERO
    scale  = ZERO
    IF (i > m) GO TO 210
  
    DO  k = i, m
       scale = scale + ABS(u(k,i))
    END DO
  
    IF (scale == ZERO) GO TO 210
  
    DO  k = i, m
       u(k,i) = u(k,i) / scale
       s = s + u(k,i)**2
    END DO
  
    f = u(i,i)
    g = -SIGN(SQRT(s),f)
    h = f * g - s
    u(i,i) = f - g
  
    DO  j = i+1, n
       s = ZERO   
       DO  k = i, m
          s = s + u(k,i) * u(k,j)
       END DO
       f = s / h
    
       DO  k = i, m
          u(k,j) = u(k,j) + f * u(k,i)
       END DO
    END DO
    u(i:m,i) = scale*u(i:m,i)
  
210 w(i)  = scale * g
    g     = ZERO
    s     = ZERO
    scale = ZERO
    IF (i > m .OR. i == n) GO TO 290
  
    DO  k = i+1, n
       scale = scale + ABS(u(i,k))
    END DO
  
    IF (scale == ZERO) GO TO 290
  
    DO  k = i+1, n
       u(i,k) = u(i,k) / scale
       s = s + u(i,k)**2
    END DO
  
    f        = u(i,i+1)
    g        = -SIGN(SQRT(s),f)
    h        = f * g - s
    u(i,i+1) = f - g
  
    DO  k = i+1, n
       rv1(k) = u(i,k) / h
    END DO
  
    IF (i == m) GO TO 270
  
    DO  j = i+1, m
       s = ZERO
       DO  k = i+1, n
          s = s + u(j,k) * u(i,k)
       END DO
       DO  k = i+1, n
          u(j,k) = u(j,k) + s * rv1(k)
       END DO
    END DO
  
270 DO  k = i+1, n
       u(i,k) = scale * u(i,k)
    END DO
  
290 anorm = MAX(anorm,ABS(w(i))+ABS(rv1(i)))
 END DO
 
!Part 3:
 IF (.NOT. matv) GO TO 410
 DO i=n,1,-1
    IF (i == n) GO TO 390
    IF (g == ZERO) GO TO 360
  
    DO  j = i+1, n
       v(j,i) = (u(i,j) / u(i,i+1)) / g
    END DO
  
    DO  j = i+1, n
       s = ZERO
       DO  k = i+1, n
          s = s + u(i,k) * v(k,j)
       END DO
       DO  k = i+1, n
          v(k,j) = v(k,j) + s * v(k,i)
       END DO
    END DO
  
360 v(i,i+1:n) = ZERO
    v(i+1:n,i) = ZERO
  
390 v(i,i) = ONE
    g      = rv1(i)
 END DO
  
!Part 4:
410 IF (.NOT. matu) GO TO 510    
  DO i=MIN(m,n),1,-1
     g = w(i)
     IF (i == n) GO TO 430
  
     DO  j = i+1, n
        u(i,j) = ZERO
     END DO
  
430  IF (g == ZERO) THEN
      u(i:m,i) = ZERO
     ELSE  
      IF (i == MIN(m,n)) GO TO 460
      DO  j = i+1, n
         s = ZERO
         DO  k = i+1, m
            s = s + u(k,i) * u(k,j)
         END DO
         f = (s / u(i,i)) / g
    
         DO  k = i, m
            u(k,j) = u(k,j) + f * u(k,i)
         END DO
      END DO
  
460   DO  j = i, m
         u(j,i) = u(j,i) / g
      END DO  
     END IF  
     
     u(i,i) = u(i,i) + ONE
  END DO
  
!Part 5:
510 DO k=n,1,-1
       k1  = k-1
       its = 0
520    DO L=k,1,-1
          IF (ABS(rv1(L)) + anorm <= anorm) GO TO 565
          IF (ABS(w(L-1)) + anorm <= anorm) EXIT
       END DO
       c = ZERO
       s = ONE
  
       DO  i = L, k
          f = s * rv1(i)
          rv1(i) = c * rv1(i)
          IF (ABS(f) + anorm <= anorm) EXIT
          g    = w(i)
          h    = SQRT(f*f+g*g)
          w(i) = h
          c    = g / h
          s    = -f / h
          IF (matu) THEN
           DO  j = 1, m
              y = u(j,L-1)
              z = u(j,i)
              u(j,L-1) = y * c + z * s
              u(j,i)   = -y * s + z * c
           END DO
          END IF
       END DO
    
!Part 6:
565 z = w(k)
    IF (L == k) GO TO 650
    IF (its == 30) GO TO 1000
    its = its + 1
    x   = w(L)
    y   = w(k1)
    g   = rv1(k1)
    h   = rv1(k)
    f   = ((y - z) * (y + z) + (g - h) * (g + h)) / (TWO * h * y)
    g   = SQRT(f*f+ONE)
    f   = ((x - z) * (x + z) + h * (y / (f + SIGN(g,f)) - h)) / x
    c   = ONE
    s   = ONE
  
    DO  i1 = L, k1
       i = i1 + 1
       g = rv1(i)
       y = w(i)
       h = s * g
       g = c * g
       z = SQRT(f*f+h*h)
       rv1(i1) = z
       c = f / z
       s = h / z
       f = x * c + g * s
       g = -x * s + g * c
       h = y * s
       y = y * c
       IF (matv) THEN
        DO  j = 1, n
           x = v(j,i1)
           z = v(j,i)
           v(j,i1) = x * c + z * s
           v(j,i)  = -x * s + z * c
        END DO
       END IF    
       z = SQRT(f*f+h*h)
       w(i1) = z

       IF (z /= ZERO) THEN
        c = f / z
        s = h / z
       END IF  
       f = c * g + s * y
       x = -s * g + c * y
       IF (matu) THEN
        DO  j = 1, m
           y = u(j,i1)
           z = u(j,i)
           u(j,i1) = y * c + z * s
           u(j,i)  = -y * s + z * c
        END DO
       END IF
    END DO
  
    rv1(L) = ZERO
    rv1(k) = f
    w(k)   = x
    GO TO 520
    
!Part 7:
650 IF (z >= ZERO) CYCLE
    w(k) = -z
    IF (matv) v(1:n,k) = -v(1:n,k)
  
 END DO
 RETURN
  
!Part 8:
1000 ierr = k
 RETURN
!*********************************************************************************************
 End
!###########################################################################################