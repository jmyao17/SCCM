        subroutine Gaussmesh()
        use vaphfb_par
        implicit none

        integer k
        real*8 xxbet(AMP%NLEG_BET)
!+---------------------------------------------------------------------+
!   Gauss-Legendre integration points and weights                      |
!+---------------------------------------------------------------------+
!  2*pi*2*pi*2 = 8*pi**2
!  wwalp,wwbet,wwgam are defined in the above regions
!  .....................................................................

          if(.NOT. ALLOCATED(sdjmk)) &
     &    ALLOCATE(sdjmk(0:AMP%jjmax,0:AMP%jjmax*2,0:AMP%jjmax*2,0:AMP%NLEG_BET,0:1))

        if(AMP%NLEG_ALP.lt. 0 .or. AMP%NLEG_BET.lt. 0 .or. AMP%NLEG_GAM.lt. 0) &
     &  stop ' Mesh points in Euler angles are not set properly '

        if(.NOT. ALLOCATED(AMP%alpha))  ALLOCATE(AMP%alpha(1:AMP%NLEG_ALP))
        if(.NOT. ALLOCATED(AMP%beta))   ALLOCATE(AMP%beta(1:AMP%NLEG_BET))
        if(.NOT. ALLOCATED(AMP%gamma))  ALLOCATE(AMP%gamma(1:AMP%NLEG_GAM))
        if(.NOT. ALLOCATED(AMP%wwalp))  ALLOCATE(AMP%wwalp(1:AMP%NLEG_ALP))
        if(.NOT. ALLOCATED(AMP%wwbet))  ALLOCATE(AMP%wwbet(1:AMP%NLEG_BET))
        if(.NOT. ALLOCATED(AMP%wwgam))  ALLOCATE(AMP%wwgam(1:AMP%NLEG_GAM))
!      ...........Gauss-Legendre integration points and weights   
        call gauleg(0.d0,2*pi, AMP%alpha,AMP%wwalp,AMP%NLEG_ALP)
        call gauleg(-1.d0,1.d0,xxbet,AMP%wwbet,AMP%NLEG_BET)
        call gauleg(0.d0,2*pi, AMP%gamma,AMP%wwgam,AMP%NLEG_GAM)

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

!          if(.NOT. ALLOCATED(ieuler_ang)) &
!     &    ALLOCATE(ieuler_ang(1:AMP%NLEG_ALP,1:AMP%NLEG_BET,1:AMP%NLEG_GAM))
!        NNLEG_ALP=NLEG_ALP
!        NNLEG_BET=NLEG_BET
!        NNLEG_GAM=NLEG_GAM
!
!  ..........................
!       do J=0,JJmax
!          weightJ(J) = (2.d0*J+1.d0)/2.d0 !(8.d0*pi) 
!       enddo
!      .........................
        return
        end
