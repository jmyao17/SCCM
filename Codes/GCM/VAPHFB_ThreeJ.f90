  subroutine PreCalc_CG()
  USE VAPHFB_Par
  implicit none
  REAL(DP) :: CG

! ......................................
! J1,J2 are j1+1/2, j2+1/2, respectively
! M1,M2 are m1+1/2, m2+1/2, starting from -J1+1 to J1
! J12,M12 are the real values.
! .....................................
  integer j1,j2,m1,m2,J12min,J12max,J12,M12
  integer J1max,J1max2

  write(*,*) ' Preparing CG coefficients ....'
  J1max=12
  J1max2=J1max*2-1

   if(HO%jmaxp5.gt.J1max) then
      J1max =HO%jmaxp5
      J1max2=HO%jmaxp5*2-1
   endif

   if(.NOT. ALLOCATED(CG_Save)) &
 & ALLOCATE(CG_Save(1:J1max,-J1max:J1max,1:J1max,-J1max:J1max, &
 &                  0:J1max2,-J1max2:J1max2)) 

      write(*,'(a,f6.3,a)')'Main: allocated CG_Save with memory:', &
     & sizeof(CG_Save)/(1024*1024*1024.0),'G byte ..'

   CG_Save(1:J1max,-J1max:J1max,1:J1max,-J1max:J1max,0:J1max2,-J1max2:J1max2) = zero 
  do J1=1,J1max !HO%jmaxp5
  do J2=1,J1max !HO%jmaxp5
  do M1=-J1+1,J1
  do M2=-J2+1,J2
    J12min = abs(J1-J2)
    J12max = abs(J1+J2-1)
    M12    = M1+M2-1
  do J12=J12min,J12max
     CG_Save(J1,M1,J2,M2,J12,M12) = CG(J1*2-1,M1*2-1,J2*2-1,M2*2-1,2*J12,2*M12)
  enddo
  enddo
  enddo
  enddo
  enddo
  end





  FUNCTION CG(j1, m1, j2, m2, J, M)
! ....................................
! all the quantum numbers are douibled  
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

  INTERFACE

  FUNCTION FacLOG(N)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  INTEGER, INTENT(IN) :: N  
  REAL(DP) ::  FacLOG
  END FUNCTION FacLOG
  END INTERFACE
  
  INTEGER, INTENT(IN) :: j1, j2, m1, m2, J, M
  REAL(DP) :: CG

  INTEGER :: H_Min, H_Max, H
  REAL(DP) :: CG_1, CG_2

  REAL(DP) :: fact1, fact2, fact3, fact4, fact5, fact6, fact7, fact8
  REAL(DP) :: fact9, fact10

  REAL(DP) :: fact11
  
  IF (  &
       & (m1 + m2   .NE. M )        .OR. &
       & (ABS(m1)   .GT. j1)        .OR. &
       & (ABS(m2)   .GT. j2)        .OR. &
       & (ABS(M)    .GT. J )        .OR. &  
       & (MOD(j1+m1,2) .NE. 0) .OR. &
       & (MOD(j2+m2,2) .NE. 0) .OR. &  
       & (MOD(J + M,2 ).NE. 0) .OR. &
       & (j1 + J  - j2 .LT. 0)       .OR. &
       & (J  - j1 + j2 .LT. 0)       .OR. &
       & (j1 + j2 - J  .LT. 0)       .OR. &
       & (j1 .LT. 0) .OR. (j2 .LT. 0) .OR. (J .LT. 0))THEN
     CG = 0 
     RETURN
  END IF

  H_Max = MIN(J-j1+j2, J+m1+m2, J+j2+m1) / 2
  H_Min = MAX(0, m1-j1, -j1+j2+m1+m2) / 2

  IF (H_Min .GT. H_Max) STOP "Factorial: Hmin Hmax error"

  fact1 = FacLOG((j1 + J - j2) / 2) 
  fact2 = FacLOG((J - j1 + j2) / 2) 
  fact3 = FacLOG((j1 + j2 - J) / 2) 
  fact4 = FacLOG((J + m1 + m2) / 2) 
  fact5 = FacLOG((J - m1 - m2) / 2) 
  fact6 = FacLOG((J + j1 + j2) / 2 + 1) 
  fact7 = FacLOG((j1 - m1)     / 2) 
  fact8 = FacLOG((j1 + m1)     / 2) 
  fact9 = FacLOG((j2 - m2)     / 2) 
  fact10= FacLOG((j2 + m2)     / 2) 

 !  IF (   fact1 .LE. 0.0 .OR. &
 !      & fact2 .LE. 0.0 .OR. &
 !      & fact3 .LE. 0.0 .OR. &
 !      & fact4 .LE. 0.0 .OR. &
 !      & fact5 .LE. 0.0 .OR. &
 !      & fact6 .LE. 0.0 .OR. &
 !      & fact7 .LE. 0.0 .OR. &
 !      & fact8 .LE. 0.0 .OR. &
 !      & fact9 .LE. 0.0 .OR. &
 !      & fact10 .LE. 0.0 ) THEN
 !    WRITE(0,*) j1, j2, J, m1, m2
 !    WRITE(0,*) fact1, fact2, fact3, fact4, fact5
 !    WRITE(0,*) fact6, fact7, fact8, fact9, fact10
 !    STOP
  !END IF

 CG_1 = fact1 + fact2 + fact3 + fact4 + fact5 &
       & - &
       &(fact6 + fact7 + fact8 + fact9 + fact10)

 CG_2 = 0.0D0

  DO H = H_Min, H_Max

     fact1 = FacLOG((J + j2 + m1) / 2 - H) 
     fact2 = FacLOG((j1 - m1) / 2 + H) 
     fact3 = FacLOG((J - j1 + j2) / 2 - H) 
     fact4 = FacLOG((J + m1 + m2) / 2 - H) 
     fact5 = FacLOG(H) 
     fact6 = FacLOG((j1 - j2 - m1 - m2) / 2 + H) 

     fact11 = fact1 + fact2 - (fact3 + fact4 + fact5 + fact6)

     CG_2 = CG_2 + &
          & (-1)**(MOD(H+(j2+m2)/2,2)) * SQRT(DBLE(J+1)) * EXP(fact11)

  END DO

  CG = exp(0.50D0 * CG_1) * CG_2

END FUNCTION CG

  FUNCTION Wigner_3j (j1, j2, J, m1, m2, M)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  
  ! Wigner 3j-Symbol 
  
  ! j1/2 j2/2 J/2 !
  ! m1/2 m2/2 M/2 !

  INTERFACE     
     FUNCTION CG(j1, m1, j2, m2, J, M)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: j1, j2, m1, m2, J, M
       REAL(DP) :: CG
     END FUNCTION CG
  END INTERFACE

  INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M

  REAL (DP) :: Wigner_3j

  Wigner_3j = (-1)**((j1-j2-M)/2) / SQRT(DBLE(J)+1.0D0) * CG(j1,m1,j2,m2,J,-M)

  RETURN

END FUNCTION Wigner_3j


FUNCTION FacLOG(N)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

  ! return LOG(N!), not N!

  INTEGER, INTENT(IN) :: N
  
  REAL(DP) ::  FacLOG

  INTEGER :: i

  FacLOG = 1.0D0

  IF ( N .EQ. 0) THEN
     FacLOG = 0.0D0
     RETURN
  END IF
 
  IF ( N .LT. 0) THEN
     PRINT *, N
     STOP "FacLOG: Impossible to calculate N! N is negative"
  END IF
  
  DO i = 1, N
     FacLOG = FacLOG * DBLE(i)
  END DO

  FacLOG = LOG(FacLOG)
  
  RETURN

END FUNCTION FacLOG

FUNCTION DBLEFacLOG(N)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  ! returns LOG(N!!), not N!!
  INTEGER, INTENT(IN) :: N
  REAL(DP) ::  DBLEFacLOG

  INTEGER :: i

  DBLEFacLOG = 1.0D0

  IF ( N .EQ. 0) THEN
     DBLEFacLOG = 0.0D0
     RETURN
  END IF
 
  IF ( N .LT. 0) THEN
     PRINT *, N
     STOP "FacLOG: Impossible to calculate N! N is negative"
  END IF

  IF (MOD(N, 2) .EQ. 0 ) THEN ! N = even
     DO i = 2, N, 2
        DBLEFacLOG = DBLEFacLOG * DBLE(i)
     END DO
  ELSE ! N = odd
     DO i = 1, N, 2
        DBLEFacLOG = DBLEFacLOG * DBLE(i)
     END DO
  END IF

  DBLEFacLOG = LOG(DBLEFacLOG)
  
  RETURN

END FUNCTION DBLEFacLOG

FUNCTION Triangle(a,b,c)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

  ! \Delta(a/2, b/2, c/2)
  ! \Delta(a/2, b/2, c/2) = sqrt( (a/2+b/2-c/2)! (a/2-b/2+c/2)! (-a/2+b/2+c/2)! / (a/2+b/2+c/2+1)! )
  ! if a/2,b/2,c/2 do not satisfy the triangle inequality, it returns zero.

  INTERFACE
     FUNCTION FacLOG(N)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, INTENT(IN) :: N
     REAL(DP) ::  FacLOG
     END FUNCTION FacLOG
  END INTERFACE

  INTEGER, INTENT(IN) :: a,b,c

  REAL(DP) :: Triangle, temp

  REAL(DP) :: F1, F2, F3, F4

  Triangle = 0.0D0
  IF (a+b-c .LT. 0 .OR. a-b+c .LT. 0 .OR. -a+b+c .LT. 0 .OR. a+b+c+1 .LT. 0) RETURN 

  IF( MOD( a+b-c,2) .NE. 0) RETURN
  IF( MOD( a-b+c,2) .NE. 0) RETURN
  IF( MOD(-a+b+c,2) .NE. 0) RETURN
  IF( MOD( a+b+c,2) .NE. 0) RETURN

  F1 = FacLOG(( a + b - c    )/2)
  F2 = FacLOG(( a - b + c    )/2)
  F3 = FacLOG((-a + b + c    )/2)
  F4 = FacLOG(( a + b + c)/2 + 1)

  temp = F1 + F2 + F3 - F4
  Triangle = EXP(0.50D0 * temp)

!  Triangle = F1 * F2 * F3 / F4

  RETURN

END FUNCTION Triangle

!  .....................................................
   FUNCTION Wigner_6j(a,b,c,d,e,f)
!  .....................................................
!  a,b,c,d,e,f are doubled
!   ( a  b  c )
!   ( d  e  f )
!  .....................................................
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

  ! Wigner 6j Symbol

  ! a b c !
  ! d e f !
 
  INTERFACE

     FUNCTION FacLOG(N)
      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: N
       REAL(DP) ::  FacLOG
     END FUNCTION FacLOG

     FUNCTION Triangle(a,b,c)
       IMPLICIT NONE
       INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: a,b,c
       REAL(DP) :: Triangle
     END FUNCTION Triangle

  END INTERFACE

 INTEGER, INTENT(IN) :: a,b,c,d,e,f

  REAL(DP) :: Wigner_6j

  REAL(DP),DIMENSION(1:4) :: Delta  
  REAL(DP),DIMENSION(1:7) :: Fact
  REAL(DP) :: W6J1, W6J2, temp

  INTEGER :: t, t_min, t_max

  Wigner_6j = 0.0D0
  
  IF ( a < 0) RETURN
  IF ( b < 0) RETURN
  IF ( c < 0) RETURN
  IF ( d < 0) RETURN
  IF ( e < 0) RETURN
  IF ( f < 0) RETURN
  IF ( a < ABS(b-c) .OR. a > b+c ) RETURN
  IF ( a < ABS(e-f) .OR. a > e+f ) RETURN
  IF ( d < ABS(b-f) .OR. d > b+f ) RETURN
  IF ( d < ABS(e-c) .OR. d > e+c ) RETURN


  Delta(1) = DBLE(Triangle(a,b,c))
  Delta(2) = DBLE(Triangle(a,e,f))
  Delta(3) = DBLE(Triangle(d,b,f))
  Delta(4) = DBLE(Triangle(d,e,c))
  
  W6J1 = Delta(1)*Delta(2)*Delta(3)*Delta(4)

!  IF ( ABS(W6J1) .LT. 1.0D-10) THEN
!     Wigner_6j = 0.0D0
!     RETURN
!  END IF

  t_max = MIN((a+b+d+e)/2, (b+c+e+f)/2, (c+a+f+d)/2)

  t_min = MAX((a+b+c)/2, (a+e+f)/2, (d+b+f)/2, (d+e+c)/2)

!  IF ( t_min .GT. t_max) THEN
!     Wigner_6j = 0.0D0
!     RETURN
!  END IF

  W6J2 = 0.0D0

  DO t = t_min, t_max

     fact(1) = FacLOG( t + (- a - b - c)/2 )
     fact(2) = FacLOG( t + (- a - e - f)/2 )
     fact(3) = FacLOG( t + (- d - b - f)/2 )
     fact(4) = FacLOG( t + (- d - e - c)/2 )
     fact(5) = FacLOG(-t + (  a + b + d + e )/2 )
     fact(6) = FacLOG(-t + (  b + c + e + f )/2 )
     fact(7) = FacLOG(-t + (  c + a + f + d )/2 )

     temp = FacLOG(t+1) - ( fact(1) + fact(2) + fact(3) + fact(4) + fact(5) + fact(6) + fact(7))

     W6J2 = W6J2 + (-1)**(t) * EXP(temp)

!     W6J2 = W6J2 + (-1)**(t) * DBLE(FacLOG(t+1)) / &
!          & (fact(1)*fact(2)*fact(3)*fact(4)*fact(5)*fact(6)*fact(7))

  END DO

  Wigner_6j = W6J1 * W6J2

  RETURN

END FUNCTION Wigner_6j
!  .....................................................


!FUNCTION Wigner_6j2(j1,j2,j3,j4,j5,j6)
!  USE PRECISION
!  IMPLICIT NONE

!  ! calcultaes Wigner 6j Symbol
!  ! from the definition of 3j symbols
!  ! not efficient. just for check

!  ! j1/2 j2/2 j3/2 !
!  ! j4/2 j5/2 j6/2 !
 
!  INTERFACE
!     FUNCTION Wigner_3j (j1, j2, J, m1, m2, M)
!       USE PRECISION
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M
!       REAL (DP) :: Wigner_3j
!     END FUNCTION Wigner_3j
!  END INTERFACE

!  INTEGER, INTENT(IN) :: j1, j2, j3, j4, j5, j6
!  REAL(DP) :: Wigner_6j2
!  REAL(DP) :: W3j1, W3j2, W3j3, W3j4
!  INTEGER :: m1, m2, m3, m4, m5, m6, S, PHASE

!  Wigner_6j2 = 0.0D0
  
!  DO m1 = -j1, j1, 2
!  DO m2 = -j2, j2, 2
!  DO m3 = -j3, j3, 2
!  DO m4 = -j4, j4, 2
!  DO m5 = -j5, j5, 2
!  DO m6 = -j6, j6, 2

!     S = (j1+j2+j3+j4+j5+j6 - m1-m2-m3-m4-m5-m6)/2
!     PHASE = (-1)**S

!     W3j1 = Wigner_3j(j1, j2, j3, m1, m2,-m3)
!     W3j2 = Wigner_3j(j1, j5, j6,-m1, m5, m6)
!     W3j3 = Wigner_3j(j4, j5, j3, m4,-m5, m3)
!     W3j4 = Wigner_3j(j4, j2, j6,-m4,-m2,-m6)

!     Wigner_6j2 = Wigner_6j2 + PHASE * W3j1*W3j2*W3j3*W3j4
     
!  END DO
!  END DO
!  END DO
!  END DO
!  END DO
!  END DO

!  RETURN

!END FUNCTION Wigner_6j2


!FUNCTION Racah_V(j1,j2,J,m1,m2,M)
!  USE PRECISION
!  IMPLICIT NONE

!  INTERFACE

!     FUNCTION Wigner_3j (j1, j2, J, m1, m2, M)
!       USE PRECISION
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M
!       REAL (DP) :: Wigner_3j
!     END FUNCTION Wigner_3j
     
!  END INTERFACE

!  INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M

!  REAL(DP) :: Racah_V

!  Racah_V = (-1)**((-j1+j2+J)/2) * Wigner_3j(j1,j2,j1,m2,m1,m2)

!  RETURN

!END FUNCTION Racah_V

FUNCTION Racah_W(a,b,c,d,e,f)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

! W(a/2 b/2 c/2 d/2; e/2 f/2)

  INTERFACE

     FUNCTION Wigner_6j(a,b,c,d,e,f)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: a,b,c,d,e,f
       REAL(DP) :: Wigner_6j
     END FUNCTION Wigner_6j
     
  END INTERFACE

  INTEGER, INTENT(IN) :: a,b,c,d,e,f
  REAL(DP) :: Racah_W

  Racah_W = 0.0D0
  IF(MOD(a+b+c+d,2) .NE. 0) RETURN

  Racah_W = (-1)**((a+b+c+d)/2) * Wigner_6j(a,b,e,d,c,f)

  RETURN

END FUNCTION Racah_W


FUNCTION Wigner_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

! j1/2 j2/2 j3/2 !
! j4/2 j5/2 j6/2 !
! j7/2 j8/2 j9/2 !

  INTERFACE

  FUNCTION Wigner_6j(a,b,c,d,e,f)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: a,b,c,d,e,f
       REAL(DP) :: Wigner_6j
     END FUNCTION Wigner_6j

  END INTERFACE

  INTEGER, INTENT(IN) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
  REAL(DP) :: Wigner_9j

  INTEGER :: x_min, x_max, x

  Wigner_9j = 0.0D0

  IF(j1 .LT. 0) RETURN
  IF(j2 .LT. 0) RETURN
  IF(j3 .LT. 0) RETURN
  IF(j4 .LT. 0) RETURN
  IF(j5 .LT. 0) RETURN
  IF(j6 .LT. 0) RETURN
  IF(j7 .LT. 0) RETURN
  IF(j8 .LT. 0) RETURN
  IF(j9 .LT. 0) RETURN

  x_min = MAX( ABS(j1-j9), ABS(j4-j8), ABS(j2-j6))
  x_max = MIN(     j1+j9,      j4+j8,      j2+j6)

  DO x = x_min, x_max

     Wigner_9j = Wigner_9j + (-1)**(x) * (x + 1.0D0) * &
          & Wigner_6j(j1,j4,j7,j8,j9, x)  * &
          & Wigner_6j(j2,j5,j8,j4, x,j6) * &
          & Wigner_6j(j3,j6,j9, x,j1,j2)
     
  END DO
  
  RETURN

END FUNCTION Wigner_9j

FUNCTION fm(I_Ang, K_Ang)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  
  INTEGER, INTENT(IN) :: I_Ang, K_Ang
  REAL(DP) :: fm

  fm = SQRT( (I_Ang + K_Ang) * (I_Ang - K_Ang + 2) / 4.0D0 )
  
END FUNCTION fm

FUNCTION fp(I_Ang, K_Ang)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

  INTEGER, INTENT(IN) :: I_Ang, K_Ang
  REAL(DP) :: fp

  fp = SQRT( (I_Ang - K_Ang) * (I_Ang + K_Ang + 2) / 4.0D0 )
  
END FUNCTION fp

FUNCTION Wigner_smalld(j,m1,m2,beta)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
! -------------------------------
! Wigner's small d-function
! d^{j/2}_{m2/2,m1/2}(\theta)
! -------------------------------
  INTERFACE


     FUNCTION FacLOG(N)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, INTENT(IN) :: N
     REAL(DP) ::  FacLOG
     END FUNCTION FacLOG
  END INTERFACE
  
  INTEGER,  INTENT(IN) :: j, m1, m2 ! doubled
  REAL(DP), INTENT(IN) :: beta      ! in radian

  REAL(DP) :: Wigner_smalld

  REAL(DP) :: fact(1:4), prefactor, mufactor
  INTEGER :: mu_min, mu_max, mu


  Wigner_smalld = 0.0D0

  fact(1) = FacLOG((j+m1)/2)
  fact(2) = FacLOG((j-m1)/2)
  fact(3) = FacLOG((j+m2)/2)
  fact(4) = FacLOG((j-m2)/2)

  prefactor = 0.50D0 * SUM(fact(1:4))
  prefactor = EXP(prefactor)

  mu_min = MAX(0, m1-m2)  /2
  mu_max = MIN(j-m2, j+m1)/2
  
  IF (mu_min .GT. mu_max) THEN
     Wigner_smalld = (mu_min - mu_max)*1.0D10
     RETURN
  END IF

  DO mu = mu_min, mu_max
     
     fact(1) = FacLOG((j-m2)/2-mu)
     fact(2) = FacLOG((j+m1)/2-mu)
     fact(3) = FacLOG(mu+(m2-m1)/2)
     fact(4) = FacLOG(mu)

     mufactor = (-1)**mu * EXP(-SUM(fact(1:4)))

     Wigner_smalld = Wigner_smalld + mufactor * (COS(0.50D0*beta))**((2*j+m1-m2)/2-2*mu) * (-SIN(0.50D0*beta))**((m2-m1)/2+2*mu)
         
  END DO

  Wigner_smalld = Wigner_smalld * prefactor

  RETURN

END FUNCTION Wigner_smalld


FUNCTION Wigner_Dmatrix(j,m1,m2,alpha,beta,gamma)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  ! note that this returns D^\ast in the notation used in nuclear theory (Takada and Ikeda (A.16))
  INTERFACE

     FUNCTION Wigner_smalld(j,m1,m2,beta)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER,  INTENT(IN) :: j, m1, m2 ! doubled
       REAL(DP), INTENT(IN) :: beta      ! in radian       
       REAL(DP) :: Wigner_smalld
     END FUNCTION Wigner_smalld
  END INTERFACE
  
  INTEGER, INTENT(IN) :: j, m1, m2
  REAL(DP), INTENT(IN) :: alpha,beta,gamma

  COMPLEX(DP) :: Wigner_Dmatrix, IUNIT = (0.0D0,1.0D0)
 
  Wigner_Dmatrix = EXP(-IUNIT*alpha*m2/2.0D0) * Wigner_smalld(j,m1,m2,beta) * EXP(-IUNIT*gamma*m1/2.0D0)

  RETURN
  
END FUNCTION Wigner_Dmatrix


FUNCTION HOB(NN, LL, n, l, n1, l1, n2, l2, Lambda, dd)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
  
  INTERFACE

     FUNCTION FacLOG(N)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: N  
       REAL(DP) ::  FacLOG
     END FUNCTION FacLOG

     FUNCTION Triangle(a,b,c)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: a,b,c
       REAL(DP) :: Triangle
     END FUNCTION Triangle

     FUNCTION DBLEFacLOG(N)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: N  
       REAL(DP) ::  DBLEFacLOG
     END FUNCTION DBLEFacLOG

     FUNCTION Wigner_9j(j1,j2,J12,j3,j4,J34,J13,J24,J)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: j1,j2,J12,j3,j4,J34,J13,J24,J
       REAL(DP) :: Wigner_9j
     END FUNCTION Wigner_9j

     FUNCTION CG(j1, m1, j2, m2, J, M)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
       INTEGER, INTENT(IN) :: j1, j2, m1, m2, J, M
       REAL(DP) :: CG
     END FUNCTION CG

!     FUNCTION Wigner_6j2(a,b,c,d,e,f)
!       USE PRECISION
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: a,b,c,d,e,f
!       REAL(DP) :: Wigner_6j2
!     END FUNCTION Wigner_6j2

  END INTERFACE

  INTEGER,  INTENT(IN) :: NN, LL, n, l, n1, l1, n2, l2, Lambda
  REAL(DP), INTENT(IN) :: dd
  REAL(DP)             :: HOB

  REAL(DP) :: PREFACTOR, FACTOR, PHASE, PHASE1
  REAL(DP) :: term1,term2,term3,term4,term5,term6,term7,term8,term9
  INTEGER :: e, EE, e1, e2, a, b, c, d, la, lb, lc, ld
  INTEGER :: ea, eb, ec, ed

!  REAL(DP) :: Wigner_6j

! equation from nucl-th/0105009v1

! major oscillator quantum numbers
  e  = 2*n +l
  EE = 2*NN+LL
  e1 = 2*n1+l1
  e2 = 2*n2+l2

  HOB = 0.0D0
  IF (e+EE - (e1+e2) .NE. 0) RETURN
  IF(MOD(l1+l2+LL+l,2) .NE. 0) RETURN 

  PHASE = (-1)**(-(l1+l2+LL+l)/2)

  PREFACTOR = 0.50D0 * ( FacLOG(n1) + FacLOG(n2) + FacLOG(NN) + FacLOG(n) + DBLEFacLOG(2*(n1+l1)+1) &
       & + DBLEFacLOG(2*(n2+l2)+1) + DBLEFacLOG(2*(NN+LL)+1) + DBLEFacLOG(2*(n+l)+1) ) &
       & - (l1+l2+LL+l)/4.0D0*LOG(2.0D0)
  PREFACTOR = EXP(PREFACTOR)

  FACTOR = 0.0D0
  
  LOOP1: DO a =  0, MIN(e1/2,EE/2) ! approx
  LOOP2: DO la = 0, MIN(e1, EE)
     ea = 2*a + la
  LOOP3: DO b =  0, MIN(e1/2, e/2)
  LOOP4: DO lb = 0, MIN(e1,  e)
     eb = 2*b + lb
     IF(ea + eb - e1 .NE. 0) CYCLE LOOP4
     IF(ABS(Triangle(2*la,2*lb,2*l1)) .LT. 1.0D-10) CYCLE LOOP4
  LOOP5: DO c =  0, MIN(e2/2,EE/2)
  LOOP6: DO lc = 0, MIN(e2, EE)
     ec = 2*c + lc
     IF(ea + ec - EE .NE. 0) CYCLE LOOP6
     IF(ABS(Triangle(2*la,2*lc,2*LL)) .LT. 1.0D-10) CYCLE LOOP6
  LOOP7: DO d =  0, MIN(e2/2, e/2)
  LOOP8: DO ld = 0, MIN(e2,  e)
     ed = 2*d + ld
     IF(ec + ed - e2 .NE. 0) CYCLE LOOP8
     IF(eb + ed -  e .NE. 0) CYCLE LOOP8
     IF(ABS(Triangle(2*lc,2*ld,2*l2)) .LT. 1.0D-10) CYCLE LOOP8
     IF(ABS(Triangle(2*lb,2*ld,2*l )) .LT. 1.0D-10) CYCLE LOOP8

     PHASE1 = (-1)**(la+lb+lc)
     term1 = 2.0D0**((la+lb+lc+ld)/2.0D0)
     term2 = dd**((ea + ed)/2.0D0)
     term3 = (1.0D0+dd)**(-(ea + eb + ec + ed)/2.0D0)

     term4 = - FacLOG(a) - FacLOG(b) - FacLOG(c) - FacLOG(d) &
          & - DBLEFacLOG(2*(a+la)+1) - DBLEFacLOG(2*(b+lb)+1) - DBLEFacLOG(2*(c+lc)+1) - DBLEFacLOG(2*(d+ld)+1)
     term4 = EXP(term4)
     term4 = term4 * (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)

     term5 = Wigner_9j(2*la,2*lb,2*l1,2*lc,2*ld,2*l2,2*LL,2*l,2*Lambda)
     term6 = CG(2*la, 0, 2*lc, 0, 2*LL, 0)
     term7 = CG(2*lb, 0, 2*ld, 0, 2*l,  0)
     term8 = CG(2*la, 0, 2*lb, 0, 2*l1, 0)
     term9 = CG(2*lc, 0, 2*ld, 0, 2*l2, 0)

     FACTOR = FACTOR +  PHASE1*term1*term2*term3*term4*term5*term6*term7*term8*term9

  END DO LOOP8
  END DO LOOP7
  END DO LOOP6
  END DO LOOP5
  END DO LOOP4
  END DO LOOP3
  END DO LOOP2
  END DO LOOP1

  ! phase convention of NPA600
  HOB = PHASE * PREFACTOR * FACTOR

  ! phase convention of Moshinsky
!  HOB = HOB * (-1)**(LL+l+Lambda)

  RETURN

END FUNCTION HOB


      SUBROUTINE CJJ (J1,J2,J3,M1,M2,M3,C)
!     CLEBSH-GORDON COEFFICIENTS                                                                                                                                                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON IOUT,IREAD,IPRI,IPREC,CAS
      COMMON /CLOG/ FLOG(100)
      C=0.0D0
      IF (M1+M2-M3) 1000,112, 1000
  112 N=J1+J2-J3
      IF(N)1000,120,120
  120 KK1=N/2
      IF (KK1+KK1-N) 1000, 130, 1000
  130 MMX=KK1
      N=J3+J1-J2
      IF(N)1000,140,140
  140 K2=N/2
      IF(K2+K2-N)1000,150,1000
  150 N=J3+J2-J1
      IF(N)1000,160,160
  160 K3=N/2
      IF(K3+K3-N)1000,170,1000
  170 K4=(J1+J2+J3)/2+1
      N=J1+M1
      IF(N)1000,180,180
  180 K5=N/2
      IF(K5+K5-N)1000,190,1000
  190 S = FLOG(KK1+1)+FLOG(K2+1)+FLOG(K3+1)-FLOG(K4+1)+FLOG(K5+1)
      N=J1-M1
      IF (N) 1000, 210, 210
  210 KK2=N/2
      IF (KK2+KK2-N) 1000, 220, 1000
  220 IF (MMX-KK2) 240, 240, 230
  230 MMX=KK2
  240 N=J2+M2
      IF(N)1000,250,250
  250 KK3=N/2
      IF (KK3+KK3-N) 1000, 260, 1000
  260 IF (MMX-KK3) 280, 280, 270
  270 MMX=KK3
  280 N=J2-M2
      IF(N)1000,290,290
  290 K3=N/2
      IF(K3+K3-N)1000,300,1000
  300 N=J3+M3
      IF(N)1000,310,310
  310 K4=N/2
      IF(K4+K4-N)1000,320,1000
  320 N=J3-M3
      IF(N)1000,330,330
  330 K5=N/2
      IF(K5+K5-N)1000,340,1000
  340 S = (S+FLOG(KK2+1)+FLOG(KK3+1)+FLOG(K3+1)                & 
     &     +FLOG(K4+1)+FLOG(K5+1))*0.5D0
      N=J3-J2+M1
      KK4=N/2
      IF (KK4+KK4-N) 1000, 345, 1000
  345 IF(N)350,360,360
  350 MMN=-KK4
      GO TO 370
  360 MMN=0
  370 N=J3-J1-M2
      KK5=N/2
      IF (KK5+KK5-N) 1000, 375, 1000
  375 IF(N)380,400,400
  380 IF (KK5+MMN) 390, 400, 400
  390 MMN=-KK5
  400 N=MMN/2
      K1=N+N-MMN
      IF(K1)410,420,420
  410 SIGN=-1.0D0
      GO TO 500
  420 SIGN=1.0D0
  500 SUM= 0.0D0
      IF (MMX-MMN) 1000, 510, 510
  510 DO 600 MU=MMN,MMX
      K1=KK1-MU
      K2=KK2-MU
      K3=KK3-MU
      K4=KK4+MU
      K5=KK5+MU
      T=S-FLOG(MU+1)-FLOG(K1+1)-FLOG(K2+1)-FLOG(K3+1) &
     &   -FLOG(K4+1)-FLOG(K5+1)
      SUM=SUM+SIGN*EXP (T)
      SIGN=-SIGN
  600 CONTINUE
      S=J3+1
      C=SUM*SQRT (S)
 1000 RETURN
      END

