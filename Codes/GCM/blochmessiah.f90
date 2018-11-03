SUBROUTINE BlochMessiah(dim,U,V,u_can,v_can,D,C,V_new, STAT)
  USE PRECISION
  IMPLICIT NONE

  ! dim should be an even number
  
  INTEGER,                             INTENT(IN)  :: dim
  COMPLEX(DP), DIMENSION(1:dim,1:dim), INTENT(IN)  :: U, V
  REAL(DP),    DIMENSION(1:dim/2),     INTENT(OUT) :: u_can, v_can
  COMPLEX(DP), DIMENSION(1:dim,1:dim), INTENT(OUT) :: D, C, V_new

  INTEGER,                             INTENT(OUT) :: STAT

  COMPLEX(DP), DIMENSION(1:dim,1:dim) :: Ud, UUd, UUdvec, UdU, UdUvec, Ubarinv, Vbarinv
  REAL(DP), DIMENSION(1:dim) :: UUdEV, UdUEV
  REAL(DP),    DIMENSION(1:dim,1:dim)  :: Ubar, Vbar
  COMPLEX(DP), DIMENSION(1:dim,1:dim) :: U_check, V_check, Check, D2

  INTEGER :: i, j,k

  REAL(DP) :: eps = 1.0D-8

  STOP "BLOCH MESSIAH under construction"
  
  ! dimension check
  IF( MOD(dim,2) .NE. 0 ) THEN
     WRITE(0,*) "BlochMessiah: dimension not 2M. exiting"
     STAT = -2
     RETURN
  END IF
  
  U_can(1:dim/2) = 0.0D0
  V_can(1:dim/2) = 0.0D0
  D(1:dim,1:dim) = (0.0D0,0.0D0)
  D2(1:dim,1:dim) = (0.0D0,0.0D0)
  C(1:dim,1:dim) = (0.0D0,0.0D0)
  UUdEV(1:dim)   = 0.0D0

  Ubar(1:dim,1:dim) = 0.0D0
  Vbar(1:dim,1:dim) = 0.0D0
  U_check(1:dim,1:dim) = 0.0D0
  V_check(1:dim,1:dim) = 0.0D0
  
  CALL UVUnitarityCheck(dim, U, V, STAT)
  IF (STAT .NE. 0) THEN
!     WRITE(0,*) "BlochMessiah: UV unitarity bad!"
!     RETURN
  END IF
     
  Ud (1:dim,1:dim) = TRANSPOSE(CONJG(U(1:dim,1:dim)))        ! U^+
  UUd(1:dim,1:dim) = MATMUL(U(1:dim,1:dim), Ud(1:dim,1:dim)) ! U U^+ = D Ubar^2 D^+
  UdU(1:dim,1:dim) = MATMUL(Ud(1:dim,1:dim), U(1:dim,1:dim)) ! U^+ U = C^+ Ubar^2 C

  CALL HermiteCheck(dim, UUd, STAT)
  IF (STAT .NE. 0) THEN
     WRITE(0,*) "BlochMessiah: UU+ not hermite"
  END IF
  CALL HermiteCheck(dim, UdU, STAT)
  IF (STAT .NE. 0) THEN
     WRITE(0,*) "BlochMessiah: U+U not hermite"
  END IF

  CALL EIGEN_ZHEEV(dim, UUd, UUdvec, UUdEV, STAT) ! UUdvec = D !
  IF (STAT .NE. 0) THEN
     WRITE(0,*) "BlochMessiah: UU+ diagonalization failed"
  END IF
  CALL EIGEN_ZHEEV(dim, UdU, UdUvec, UdUEV, STAT)
  IF (STAT .NE. 0) THEN
     WRITE(0,*) "BlochMessiah: U+U diagonalization failed"
  END IF

  ! CHECK diagonalization is performed correctly
  CHECK(1:dim,1:dim) = 0.0D0
  DO i = 1, dim
  DO j = 1, dim
  DO k = 1, dim
     CHECK(i,j) = CHECK(i,j) + UUdvec(i,k)*UUdEV(k)*CONJG(UUdvec(j,k))  !+ CONJG(UUdvec(k,i))*UUdEV(k)*UUdvec(k,j)
  END DO
  END DO
  END DO

  DO i = 1, dim
  DO j = 1, dim
     IF( ABS(CHECK(i,j) - UUd(i,j)) .GT. 1.0D-8) THEN
        PRINT *, i, j, UUd(i,j), CHECK(i,j)
     END IF
  END DO
  END DO
!  STOP "DEBUG"
  DO i = 1, dim/2
     ! All eigenvalue must be paired
     IF(ABS(UUdEV(2*i-1)-UUdEV(2*i)) .GT. 1.0D-10) THEN
        WRITE(0,*) "BlochMessiah: Ubar from UU+ wrong"
        WRITE(0,*) "u_i", UUdEV(2*i-1), UUdEV(2*i)
        STAT = -10
        RETURN
     END IF

     IF(ABS(UdUEV(2*i-1)-UdUEV(2*i)) .GT. 1.0D-10) THEN
        WRITE(0,*) "BlochMessiah: Ubar from U+U wrong"
        WRITE(0,*) "u_i", UdUEV(2*i-1), UdUEV(2*i)
        STAT = -20
        RETURN
     END IF

     ! And the eigenvalue from UU+ and U+U should be the same
     IF(ABS(UUdEV(2*i-1)-UdUEV(2*i-1)) .GT. 1.0D-10) THEN
        WRITE(0,*) "BlochMessiah: Ubar from U+U and UU+ do not agree"
        WRITE(0,*) "u_i", UUdEV(2*i-1), UdUEV(2*i-1)
        STAT = -25
        RETURN
     END IF
     
     ! eigenvalue range u_i =< 1
     IF( UUdEV(2*i-1) .GT. 1.0D0+1.0D-8) THEN
        WRITE(0,*) "BlochMessiah: u_i from UU+ larger than 1"
        WRITE(0,*) "u_i", UUdEV(2*i-1)
        STAT = -30
        RETURN
     END IF
     IF( UdUEV(2*i-1) .GT. 1.0D0+1.0D-8) THEN
        WRITE(0,*) "BlochMessiah: u_i from U+U larger than 1"
        WRITE(0,*) "u_i", UdUEV(2*i-1)
        STAT = -40
        RETURN
     END IF
  
     ! remove numerical error
     IF( UUdEV(2*i-1) .GT. 1.0D0 - eps .OR. UUdEV(2*i) .GT. 1.0D0 - eps) THEN
        UUdEV(2*i-1) = 1.0D0 - eps**2
        UUdEV(2*i)   = 1.0D0 - eps**2
     END IF
     IF( UdUEV(2*i-1) .GT. 1.0D0 - eps .OR. UdUEV(2*i) .GT. 1.0D0 - eps ) THEN
        UdUEV(2*i-1) = 1.0D0 - eps**2
        UdUEV(2*i)   = 1.0D0 - eps**2
     END IF

     ! eigenvalue range u_i >= 0 
     IF( UUdEV(2*i-1) .LT. -1.0D-8) THEN
        WRITE(0,*) "BlochMessiah: u_i from UU+ smaller than 0"
        WRITE(0,*) "u_i", UUdEV(2*i-1)
        STAT = -50
        RETURN
     END IF
     IF( UdUEV(2*i-1) .LT. -1.0D-8) THEN
        WRITE(0,*) "BlochMessiah: u_i from U+U smaller than 0"
        WRITE(0,*) "u_i", UdUEV(2*i-1)
        STAT = -60
        RETURN
     END IF
  
     ! remove numerical error
     IF( UUdEV(2*i-1) .LT. eps .OR. UUdEV(2*i) .LT. eps ) THEN
        UUdEV(2*i-1) = eps**2
        UUdEV(2*i)   = eps**2
     END IF
     IF( UdUEV(2*i-1) .LT. eps .OR. UdUEV(2*i) .LT. eps ) THEN
        UdUEV(2*i-1) = eps**2
        UdUEV(2*i)   = eps**2
     END IF

     u_can(i) = SQRT(UUdEV(2*i))
     v_can(i) = SQRT(1.0D0 - UUdEV(2*i))
     PRINT *, "i = ", i, u_can(i), v_can(i)
     

     Ubar(2*i-1,2*i-1) =  SQRT(UUdEV(2*i))
     Ubar(2*i,2*i)     =  SQRT(UUdEV(2*i))
     Vbar(2*i-1,2*i)   =  SQRT(1.0D0-UUdEV(2*i))
     Vbar(2*i,2*i-1)   = -SQRT(1.0D0-UUdEV(2*i))

  END DO

  PRINT *, "summation over the states"
  PRINT *, " u_i^2      : ", SUM(u_can(1:dim/2)**2)
  PRINT *, " v_i^2      : ", SUM(v_can(1:dim/2)**2)
  PRINT *, " u_i^2+v_i^2: ", SUM(u_can(1:dim/2)**2+v_can(1:dim/2)**2)
  
  D2(1:dim,1:dim) = UUdvec(1:dim,1:dim)
  !  C(1:dim,1:dim) = UdUvec(1:dim,1:dim)

  Ubarinv(1:dim,1:dim) = 0.0D0
  DO i = 1, dim
     PRINT *, "u_can ", i, Ubar(i,i), Vbar(i,i+1)
     Ubarinv(i,i) = Ubar(i,i)**(-1)
  END DO

  C(1:dim,1:dim) = MATMUL( Ubarinv(1:dim,1:dim), MATMUL( CONJG(TRANSPOSE(D2(1:dim,1:dim))), U(1:dim,1:dim)))

  ! Unitarity check
  Check(1:dim,1:dim) = MATMUL( TRANSPOSE(CONJG(C(1:dim,1:dim))), C(1:dim,1:dim))
  PRINT *, "Unitarity check of C"
  DO i = 1, dim
  DO j = 1, dim
     IF( i .EQ. j .AND. ABS(Check(i,j)-1.0D0) .GT. 1.0D-7) THEN
        PRINT *, i, j, Check(i,j)
     END IF
     IF( i .NE. j .AND. ABS(Check(i,j)) .GT. 1.0D-7) THEN
        PRINT *, i, j, Check(i,j)
     END IF
  END DO
  END DO
  
  Check(1:dim,1:dim) = MATMUL( TRANSPOSE(CONJG(D2(1:dim,1:dim))), D2(1:dim,1:dim))
  PRINT *, "Unitarity check of D"
  DO i = 1, dim
  DO j = 1, dim
     IF( i .EQ. j .AND. ABS(Check(i,j)-1.0D0) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
     IF( i .NE. j .AND. ABS(Check(i,j)) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
  END DO
  END DO

!  Vbarinv(1:dim,1:dim) = 0.0D0
!  DO i = 1, dim/2     
!     Vbarinv(2*i-1,2*i) =  Vbar(2*i,2*i-1)**(-1)
!     Vbarinv(2*i,2*i-1) =  Vbar(2*i-1,2*i)**(-1)
!  END DO

  ! calc V part
  
  Vbar(1:dim,1:dim) = 0.0D0

  Vbar(1:dim,1:dim) = MATMUL( TRANSPOSE(D2(1:dim,1:dim)), MATMUL( V(1:dim,1:dim), CONJG(TRANSPOSE(C(1:dim,1:dim)))))
  
  !CHECK Vbarinv
  PRINT *, "Vbar from D and C of U"

  DO i = 1, dim
  DO j = 1, dim
     IF( ABS(Vbar(i,j)) .GT. 1.0D-10) PRINT *, i, j, Vbar(i,j)        
  END DO
END DO
  
  
  Check(1:dim,1:dim) = MATMUL( Vbar(1:dim,1:dim), Vbarinv(1:dim,1:dim))
  PRINT *, "Orthogonality check of Vbar Vbar^{-1}"
  DO i = 1, dim
  DO j = 1, dim
     IF( i .EQ. j .AND. ABS(Check(i,j)-1.0D0) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
     IF( i .NE. j .AND. ABS(Check(i,j)) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
  END DO
  END DO  
  
!  D(1:dim,1:dim) = CONJG( MATMUL( V(1:dim,1:dim), &
!       & MATMUL( CONJG(TRANSPOSE(C(1:dim,1:dim))), Vbarinv(1:dim,1:dim))))
  D = D2
  Check(1:dim,1:dim) = MATMUL( TRANSPOSE(CONJG(D(1:dim,1:dim))), D(1:dim,1:dim))
  PRINT *, "Unitarity check of D new"
  DO i = 1, dim
  DO j = 1, dim
     IF( i .EQ. j .AND. ABS(Check(i,j)-1.0D0) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
     IF( i .NE. j .AND. ABS(Check(i,j)) .GT. 1.0D-8) THEN
        PRINT *, i, j, Check(i,j)
     END IF
  END DO
  END DO

!  PRINT *, "Dmat"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (ABS(D(i,j)), j = 1, 8)
!  END DO
!  PRINT *, "D2mat"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (ABS(D2(i,j)), j = 1, 8)
!  END DO
  !  PRINT *, "Cmat"
!  DO i = 1, 4
!     WRITE(*,'(8ES15.5)') (C(i,j), j = 1, 4)
!  END DO
!  PRINT *, "Ubar"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (Ubar(i,j), j = 1, 8)
!  END DO
!  PRINT *, "Vbar"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (Vbar(i,j), j = 1, 8)
!  END DO
!  PRINT *, "U"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (U(i,j), j = 1, 8)
!  END DO
!  PRINT *, "V"
!  DO i = 1, 8
!     WRITE(*,'(8ES15.5)') (V(i,j), j = 1, 8)
!  END DO 
  
  ! check U and V from the Bloch-Messiah decomposition

  U_check(1:dim,1:dim) = MATMUL(       D(1:dim,1:dim),  MATMUL( Ubar(1:dim,1:dim), C(1:dim,1:dim)) )
  V_check(1:dim,1:dim) = MATMUL( CONJG(D(1:dim,1:dim)), MATMUL( Vbar(1:dim,1:dim), C(1:dim,1:dim)) )
  V_new(1:dim,1:dim) = V_check(1:dim,1:dim)

  PRINT *, "CHECKING U matrix from Bloch-Messiah"
  DO i = 1, dim
  DO j = 1, dim
     IF( ABS(U_check(i,j)  - U(i,j)) .GT. 1.0D-10) THEN
        PRINT *, i, j, U(i,j), U_check(i,j)
     END IF
  END DO
  END DO
  PRINT *, "CHECKING V matrix from Bloch-Messiah -always failed"
  DO i = 1, dim
  DO j = 1, dim
     IF( ABS(V_check(i,j)  - V(i,j)) .GT. 1.0D-8) THEN
        PRINT *, i, j, V(i,j), V_check(i,j)
     END IF
  END DO
  END DO

!  PRINT *, "UV unitarity check1"
!  CALL UVUnitarityCheck(dim, U, V, STAT)
!  PRINT *, "UV unitarity check2"
!  CALL UVUnitarityCheck(dim, U, V_check, STAT)
!  PRINT *, "UV unitarity check3"
!  CALL UVUnitarityCheck(dim, U_check, V, STAT)
  PRINT *, "UV unitarity check4"
  CALL UVUnitarityCheck(dim, U_check, V_check, STAT)
  PRINT *, "Check done: STAT = ", STAT
  
  PRINT *, "END SUBROUTINE Bloch-Messiah"
  RETURN

END SUBROUTINE BlochMessiah
  
