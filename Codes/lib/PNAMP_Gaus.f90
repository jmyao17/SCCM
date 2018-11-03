        subroutine Gaussmesh()
!+---------------------------------------------------------------------+
!   Gauss-Legendre integration points and weights                      |
!   2*pi*2*pi*2 = 8*pi**2
!   wwalp,wwbet,wwgam are defined in the above regions
!+---------------------------------------------------------------------+
        use vaphfb_par
        implicit none

        integer k
        real*8 xxbet(AMP%NLEG_BET)
!  .....................................................................

        if(AMP%NLEG_ALP.lt. 0 .or. AMP%NLEG_BET.lt. 0 .or. AMP%NLEG_GAM.lt. 0) &
     &  stop ' Mesh points in Euler angles are not set properly '

        if(.NOT. ALLOCATED(AMP%alpha))  ALLOCATE(AMP%alpha(1:AMP%NLEG_ALP))
        if(.NOT. ALLOCATED(AMP%beta))   ALLOCATE(AMP%beta(1:AMP%NLEG_BET))
        if(.NOT. ALLOCATED(AMP%gamma))  ALLOCATE(AMP%gamma(1:AMP%NLEG_GAM))
        if(.NOT. ALLOCATED(AMP%wwalp))  ALLOCATE(AMP%wwalp(1:AMP%NLEG_ALP))
        if(.NOT. ALLOCATED(AMP%wwbet))  ALLOCATE(AMP%wwbet(1:AMP%NLEG_BET))
        if(.NOT. ALLOCATED(AMP%wwgam))  ALLOCATE(AMP%wwgam(1:AMP%NLEG_GAM))
!      ...........Gauss-Legendre integration points and weights   
       if(AMP%i3DAMP.eq.0) then
          ! 1DAMP
          !call gauleg(0.d0,2*pi, AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
          !call gauleg(0.d0,2*pi, AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
          call gauleg(-1.d0,1.d0,xxbet,AMP%wwbet,AMP%NLEG_BET)
          call trapezoidal(0.d0,2*pi,AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
          call trapezoidal(0.d0,2*pi,AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
        else if(AMP%i3DAMP.ne.0) then
          ! 3DAMP for triaxial without cranking
          if(Input%icr.eq.0) then
          ! 3DAMP without cranking
             !call gauleg(0.d0,pi/2.0, AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
             !call gauleg(0.d0,pi, AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
             call gauleg(0.d0,1.d0,xxbet,AMP%wwbet,AMP%NLEG_BET)
             call trapezoidal(0.d0,pi/2.0,AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
             call trapezoidal(0.d0,pi,AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
          else if(Input%icr.eq.1) then
          ! 3DAMP with cranking
             call gauleg(0.d0,1.d0,xxbet,AMP%wwbet,AMP%NLEG_BET)
             call trapezoidal(0.d0,pi,AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
             call trapezoidal(0.d0,pi*2,AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
           !  call gauleg(0.d0,pi*2, AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
          !   call gauleg(0.d0,pi, AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
           else  ! Input%icr.eq.2
          ! .................................
          ! general case(no symmetry)
          ! alpha,gamma in between [0, 2pi]
          ! beta        in between [0, pi]
          ! .................................
          ! call gauleg(0.d0,pi*2, AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
          !  call gauleg(0.d0,pi*2, AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)

            call gauleg(-1.d0,1.d0,xxbet,AMP%wwbet,AMP%NLEG_BET)
            call trapezoidal(0.d0,2*pi,AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)
            call trapezoidal(0.d0,2*pi,AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
           endif
        endif
        do k=1,AMP%NLEG_BET
           AMP%beta(k) = dacos(xxbet(k))
        end do
!  ......................... a special case
        if(AMP%NLEG_ALP .eq. 1) then
           AMP%alpha(1)=0.d0
           AMP%wwalp(1) =1.d0
        endif
        if(AMP%NLEG_BET.eq. 1) then
           AMP%beta(1) =0.d0
           AMP%wwbet(1) =1.d0
        endif
        if(AMP%NLEG_GAM .eq. 1) then
           AMP%gamma(1)=0.d0
           AMP%wwgam(1) =1.d0
        endif
!      .........................
        return
        end


      subroutine trapezoidal(x1,x2,x,w,n)
      INTEGER n,m
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION xL  

        xL = x2 -x1
        DO m=1,n                   
          x(m)  = (m-0.5)*xL/n            
          w(m)  = 1.d0/n  
        enddo

      return
      END

