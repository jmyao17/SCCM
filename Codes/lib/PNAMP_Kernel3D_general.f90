        subroutine Kernel4general(lpr,Nphi1,Nphi2,NLEV)
        USE VAPHFB_PAR
        ! .................................
        ! computing <q0| O P^J_{k,k'}| q1> 
        ! .................................
        ! Note: the Varshalovich notation for the rotation operator:
        ! R(alp,bet,gam)   = e^(-i alp J_z) * e^(-i bet J_y) * e^(-i gam J_z)
        ! D^J*_(K1,K2)     = e^(i alp K1) * d^J_(K1,K2) (bet) * e^(i gam K2)
        ! ...................................................................
        ! ATTENTION:
        ! The small d-function: d^J_(K1,K2) (bet) is different from that
        ! of Edmond's notation in the following way:
        ! d^J_(K1,K2) (bet) here = d^J_(K1,K2) (-bet) by Edmonds
        ! ...................................................................
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)
        logical lpr
        Integer a12,a34,b12,Ncosbet                                     ! symmetry for beta angle 
        Integer ieuler_ang(1:AMP%NLEG_ALP,1:AMP%NLEG_BET,1:AMP%NLEG_GAM)
        real*8  start1,finish1,start,finish

        !WAVE FUNCTIONS
        DIMENSION UAUX(NLEV,NLEV)  !auxiliary vacuum
        DIMENSION VAUX(NLEV,NLEV)  !auxiliary vacuum
        DIMENSION ZU0(NLEV,NLEV) !q (bra)
        DIMENSION ZV0(NLEV,NLEV) !q (bra)
        DIMENSION ZU1(NLEV,NLEV) !q' (ket)
        DIMENSION ZV1(NLEV,NLEV) !q' (ket)
        DIMENSION ZUAUX(NLEV,NLEV) !q (bra)
        DIMENSION ZVAUX(NLEV,NLEV) !q (bra)
        DIMENSION ZU1bar(NLEV,NLEV) !q' (ket)
        DIMENSION ZV1bar(NLEV,NLEV) !q' (ket)
        Dimension ZZTZ(NLEV,NLEV)
        dimension  hjkk(JJmax,-JJmax:JJmax,-JJmax:JJmax)
        !DENSITIES AND FIELDS
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION zkapa10(NLEV,NLEV)       ! pairing tensor 10 (cc)
        DIMENSION zkapa01(NLEV,NLEV)       ! pairing tensor 01 (c+c+)
        DIMENSION zgamma(NLEV,NLEV) ! HF field
        DIMENSION zham(NLEV,NLEV)   ! kinetic+HF field
        DIMENSION zdelta10(NLEV,NLEV)       ! pairing field 10
        DIMENSION zdelta01(NLEV,NLEV)       ! pairing field 10
!      ................ 1B,2B,densities in J-scheme
        DIMENSION ZOBDJ0    (0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        DIMENSION ZRho1BJ0  (0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        DIMENSION ZLambda2B (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        Complex*16, allocatable, DIMENSION(:,:) :: ZRO3BJJ,ZRho3BJJ,ZLambda3B    
        real*8, allocatable, DIMENSION(:,:) :: aa    
        real*8, Allocatable, DIMENSION(:) :: EJ,FK    

        REAL*8, DIMENSION(0:1,1:HO%nMax,0:HO%lmax*2,0:1,1:HO%nMax,0:HO%LMax*2) :: Rho1B_Full 


        !MATRIX ELEMENTS
        DIMENSION ZROT(NLEV,NLEV)
        DIMENSION ZAM_phi_omega(NLEV,NLEV)

        !OTHER QUANTITIES
        DIMENSION ZZc0(NLEV,NLEV),ZZc1(NLEV,NLEV)
        DIMENSION Zsumphi(NOBS)
        DIMENSION ZRES(AMP%NNNP,NOBS)
        DIMENSION ZsumJKK(NOBS,-JJmax-2:JJmax+2,-JJmax:JJmax)
        CHARACTER*100 cwf,cwf1,cwf2
        complex*16 zqlm(-2:2,2),zqlm2(1:4,-2:2,2)
        complex*16 zqlmc(-2:2,2)
        complex*16 ztest,zqJ0(-2:2,2)

        real*8 t0,t1,t2,t3,t4,t5,t6,t7,t8 !,omp_get_wtime()

      print *, '3j'
      print *, wignei(2,2,0,2,-2,0)
      print *, wignei(2,2,0,0,0,0)

        t0 = omp_get_wtime()
        write(*,*) ' -> Start to compute kernels for arbit. states ..'
!  ....................... preparation
        Aneut = Input%nneut !Nucl%nucleon(0) - Nucl%ncore(0)
        Aprot = Input%nprot !Nucl%nucleon(1)
        aaa   = 1.d0*NLEV*(NLEV+1)/2.d0
        snpfaf= (-1.d0)**(aaa)         !for the Pfaffian, check it!
        NLEG_BET_Max = AMP%NLEG_BET

        Ncosbet = 1  
!  ..........................
        do J=0,JJmax
           !weightJ(J) = (2.d0*J+1.d0)/(8.d0*pi**2)  
           weightJ(J) = (2.d0*J+1.d0)/2.d0  
        enddo
      !read wave functions
      if(.NOT. ALLOCATED(HFB%U0))  ALLOCATE(HFB%U0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%V0))  ALLOCATE(HFB%V0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%U1))  ALLOCATE(HFB%U1(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%V1))  ALLOCATE(HFB%V1(1:HO%NLEV,1:HO%NLEV))
      !read u0,v0 of ref. state 1
          cwf1=File%cwf(Nphi1)
          print *, cwf1
          open(unit=10,file=cwf1)
          read(10,*) HFB%ide0(0),HFB%ide0(1)
          do ii=1,NLEV
           do jj=1,NLEV
            read (10,*)  HFB%U0(ii,jj),HFB%V0(ii,jj)
           end do
          end do
          close(10)
          do ii=1,NLEV
           do jj=1,NLEV
            ZU0(ii,jj)=HFB%U0(ii,jj)*zone
            ZV0(ii,jj)=HFB%V0(ii,jj)*zone
          end do
          end do


      !read u1,v1 of ref. state 2
          cwf2=File%cwf(Nphi2)
          print *, cwf2
          open(unit=11,file=cwf2)
          read(11,*) HFB%ide1(0),HFB%ide1(1)
          do ii=1,NLEV
           do jj=1,NLEV
              read (11,*)  HFB%U1(ii,jj),HFB%V1(ii,jj)
           end do
          end do
          close(11)
          do ii=1,NLEV
           do jj=1,NLEV
            ZU1(ii,jj)=HFB%U1(ii,jj)*zone
            ZV1(ii,jj)=HFB%V1(ii,jj)*zone
           end do
          end do
!     check if anyone of the two wavefunctions is HF state
!     ide(0:1) = 0 (HF); 1 (HFB)
      do it=0,1
         HFB%ide(it)=HFB%ide0(it)*HFB%ide1(it)
      enddo !it

       
      ! check if densities are to be calculated 
            open(69,file=File%Rho1B,status='unknown')
         if(Input%idens.ne.0) then 
            open(300,file=File%Rho2B,status='unknown')
            if(Input%iRho3B.eq.1) &
            open(89,file=File%Rho3B,status='unknown')
         endif


!        ......... for shell-model space, always using PNP
          if(Input%IntType.eq.0) then
            HFB%ide(0) = 1 
            HFB%ide(1) = 1
         endif
!        ........ if there is pairing in both neutrons and protons, use
!        pfaffian method 
!        ......... 10^N  * 10^N is from pfaffian
         if(HFB%ide(0)+HFB%ide(1) .ge.1) then  
            call zbareover(ZV0,ZU0,ZZc0,Det0,NLEV)  !q (bra)
            zDetvac0=Det0*10.d0**(-NLEV*Input%iscale/2) !a bug is removed here, zone*((1./Det0)**.25)
            call zbareover(ZV1,ZU1,ZZc1,Det1,NLEV)  !q' (ket)
            zDetvac1=Det1*10.d0**(-NLEV*Input%iscale/2) !zone*((1./Det1)**.25)
            !write(*,*) '<phi|phi>(0):',zDetvac0,zDetvac1

            PNP%NFOM = PNP%MPhi 
          else
            PNP%NFOM = 1
            print *, 'Caution: At least, one of the states is HF state!'
            print *, '         Onishi formula is used for the norm!'
         endif
!    ........................... initialization
        do ii=1,NOBS
         do jj=1,AMP%NNNP
            ZRES(jj,ii)=zzero
         end do
        end do
!    ......... density
      ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)       = zzero

        if(Input%idens.ne.0) then
           if(.NOT. ALLOCATED(Dens%ZRho2BJJ))  &
           allocate(Dens%ZRho2BJJ(0:aMaxMax,0:aMaxMax,0:bMax))
           if(.NOT. ALLOCATED(Dens%ZTBDJJ))  &
           allocate(Dens%ZTBDJJ(0:aMaxMax,0:aMaxMax,0:bMax))
           if(Input%iRho3B.eq.1) then
               if(.NOT. ALLOCATED(ZRho3BJJ))  &
     &         ALLOCATE(ZRho3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))

               if(.NOT. ALLOCATED(ZRO3BJJ))  &
     &         ALLOCATE(ZRO3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))

               ZRho3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax)  = zzero
               print *, 'Size of ZRho3BJJ:', &
     &         sizeof(ZRho3BJJ)/(1024*1024*1024.0),'G byte ..'
           endif

           print *, 'Size of ZRho2BJJ:', &
     &     sizeof(Dens%ZRho2BJJ)/(1024*1024*1024.0),'G byte ..'
     
            Dens%ZRho2BJJ(:,:,:) = zzero
       endif
!  .......................................
        if(.NOT. ALLOCATED(AMP%ZROT_m1m2))  ALLOCATE(AMP%ZROT_m1m2(1:HO%NLEV,1:HO%NLEV))
        print *, ' --> start loops over rotation angles'
        icheck = AMP%NLEG_BET*AMP%NLEG_ALP*AMP%NLEG_GAM/100        
        iang=0
! .......................... loop over beta angle
        DO k=1,max(NLEG_BET_Max,1)
           bet=AMP%beta(k)
           zwb=zone*AMP%wwbet(k)

           if(AMP%NLEG_BET.eq.1) then
              bet = zero
              zwb = zone
           endif
! .......................... loop over alpha angle
        DO kkk=1,AMP%NLEG_ALP
           alp=AMP%alpha(kkk)
           zwa=AMP%wwalp(kkk)
! .......................... loop over gamma angle
        DO kk=1,AMP%NLEG_GAM
           gam=AMP%gamma(kk)
           zwg=AMP%wwgam(kk)

           iang=iang+1
           ! check the progress of the calculation
           if(mod(iang,icheck).eq.0 .and. lpr) call progress_bar(int(iang/icheck)) 

           ! rotation angles
           ieuler_ang(kkk,k,kk) = iang

           ! weights
           zffJ0=zwa*zwb*zwg !*Ncosbet ! only for the DMEs of g.s. (J=0) 

!       ........... TWO-BODY Density coupled to J=0
!        The two-body density requires large memory, may cause
!        fragmentation error

          !print *, ' Initializing one-body density for a given rotation angle ...'
          ZOBDJ0 (0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)   = zzero
          if(Input%idens.ne.0) then
            !t0 = omp_get_wtime()
            !print *, ' Initializing two-body density ...'
            Dens%ZTBDJJ (:,:,:) = zzero          
            if(Input%iRho3B.eq.1)  then 
              !print *, ' Initializing three-body density ...'
              ZRO3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax)  = zzero
            endif
            !t1 = omp_get_wtime()
            !print *, ' Time consummed for initializing densities:',t1-t0
          endif
       ! Rotation matrix  
       call ROTATION_MAT(alp,bet,gam,ZROT,NLEV)
       do ia=1,NLEV
       do ib=1,NLEV
          AMP%ZROT_m1m2(ia,ib) = ZROT(ia,ib)
       enddo
       enddo
       ! initializing the variables that are used to store the
       ! number-projected quantities
       do ii=1,NOBS  !AQUI
          zsumphi(ii)=Zzero
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            FOMENKO STARTS                      !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        t0 = omp_get_wtime()
        !print *, 'Number of mesh points for PNP:', PNP%NFOM
        DO LLLP=0,PNP%NFOM-1   !Fomenko                

          phiLP=(2*pi/PNP%NFOM)*LLLP             ! phi_p
          if(PNP%NFOM.eq.1) phiLP=0.d0 

          zphas_p=cdexp(-Zimag*phiLP*Aprot)  ! exp(-i*Z*phi_p)

        DO LLLN=0,PNP%NFOM-1   !Fomenko                

           !write(*,'(a,i4,i4,a)') '  Gauge angles: (',LLLP,LLLN,')'

          phiLN=(2*pi/PNP%NFOM)*LLLN               ! phi_n
          if(PNP%NFOM.eq.1) phiLN=0.d0

          zphas_n=cdexp(-Zimag*phiLN*Aneut)     ! exp(-i*N*phi_n)

          PNP%phi(1) = phiLN
          PNP%phi(2) = phiLP
          PNP%zemiphi(1) = zphas_n
          PNP%zemiphi(2) = zphas_p
!         ...............
!         rotated U and V in gauge angles
          call UVphiomega(phiLP,phiLN,ZROT,ZU1,ZV1,ZU1bar,ZV1bar,ZAM_phi_omega,NLEV)
         !norm overlap   
         if(PNP%NFOM .eq.1) then
!           Onishi formula is not valid for cranking states?
!           sign problem: K. Enami, K. Tanabe, and N. Yoshinaga, PRC59, 135 (1999)
            call Onishi(ZU0,ZV0,ZU1bar,ZV1bar,znorm,HO%NLEV)
         else
           ! Pfaffian method
            call Pfaf_Over(ZZc0,ZU1bar,ZV1bar,ZdetPfaf,NLEV)
!            print *, snpfaf,zDetvac0*zDetvac1,ZdetPfaf,zphas_p*zphas_n
             znorm = snpfaf*zDetvac0*zDetvac1*ZdetPfaf*zphas_p*zphas_n
         endif
        zsumphi(1)=zsumphi(1)+znorm

        ! computing mixed/transition densities
        call ROKAPPAphiomega(ZAM_phi_omega,ZU0,ZV0,ZU1,ZV1,ZRO,ZKAPA10,ZKAPA01,NLEV)

        call ZN_0(ZRO,Zprot_0,Zneut_0,NLEV)
        zsumphi(2)=zsumphi(2)+znorm*Zprot_0
        zsumphi(3)=zsumphi(3)+znorm*Zneut_0


        call HFB_Energy_COMPLEX(zro,zkapa10,zkapa01,H%zME1BM,zEkin,zEHFB,NLEV,.false.)


        zsumphi(4)=zsumphi(4)+znorm*zEkin
        zsumphi(5)=zsumphi(5)+znorm*zEHF
        zsumphi(6)=zsumphi(6)+znorm*zEPa
        zsumphi(7)=zsumphi(7)+znorm*(zEHFB+H%E0)

!        write(*,*) ' Start ZJ2_0 ....'
        call ZJ2_0(ZRO,ZKAPA10,ZKAPA01,zJ2_MV,NLEV)
        zsumphi(8)=zsumphi(8)+znorm*zJ2_MV

!   ... quadrupole moments
        if(Input%iE2.eq.1) then
!        write(*,*) ' Start ZQ2 ....'
!       ...............................................
        call ZQ2(cME1B%Q22t*zone,ZRO,zQ22p,zQ22n,NLEV)
        call ZQ2(cME1B%Q21t*zone,ZRO,zQ21p,zQ21n,NLEV)
        call ZQ2(cME1B%Q20t*zone,ZRO,zQ20p,zQ20n,NLEV)
        call ZQ2(cME1B%Q2_1t*zone,ZRO,zQ2_1p,zQ2_1n,NLEV)
        call ZQ2(cME1B%Q2_2t*zone,ZRO,zQ2_2p,zQ2_2n,NLEV)

          ! .......... shell-model interactions 
          if(Input%IntType.eq.0) then
           zsumphi(9)=zsumphi(9)  +znorm*HO%b_osc**2*(Input%ep*zQ2_2p+Input%en*zQ2_2n)
           zsumphi(10)=zsumphi(10)+znorm*HO%b_osc**2*(Input%ep*zQ2_1p+Input%en*zQ2_1n)
           zsumphi(11)=zsumphi(11)+znorm*HO%b_osc**2*(Input%ep*zQ20p+Input%en*zQ20n)
           zsumphi(12)=zsumphi(12)+znorm*HO%b_osc**2*(Input%ep*zQ21p+Input%en*zQ21n)
           zsumphi(13)=zsumphi(13)+znorm*HO%b_osc**2*(Input%ep*zQ22p+Input%en*zQ22n)

          else
           zsumphi(9)=zsumphi(9)+znorm*zQ2_2p
           zsumphi(10)=zsumphi(10)+znorm*zQ2_1p
           zsumphi(11)=zsumphi(11)+znorm*zQ20p
           zsumphi(12)=zsumphi(12)+znorm*zQ21p
           zsumphi(13)=zsumphi(13)+znorm*zQ22p
          endif

!      ............ trans. of Q_2mu
        call ZQ2t(cME1B%Q22t*zone,ZRO,zQ22pt,zQ22nt,NLEV)
        call ZQ2t(cME1B%Q21t*zone,ZRO,zQ21pt,zQ21nt,NLEV)
        call ZQ2t(cME1B%Q20t*zone,ZRO,zQ20pt,zQ20nt,NLEV)
        call ZQ2t(cME1B%Q2_1t*zone,ZRO,zQ2_1pt,zQ2_1nt,NLEV)
        call ZQ2t(cME1B%Q2_2t*zone,ZRO,zQ2_2pt,zQ2_2nt,NLEV)

          if(Input%IntType.eq.0) then
             zsumphi(14)=zsumphi(14)+znorm*HO%b_osc**2*(Input%ep*zQ2_2pt+Input%en*zQ2_2nt)
             zsumphi(15)=zsumphi(15)+znorm*HO%b_osc**2*(Input%ep*zQ2_1pt+Input%en*zQ2_1nt)
             zsumphi(16)=zsumphi(16)+znorm*HO%b_osc**2*(Input%ep*zQ20pt+Input%en*zQ20nt)
             zsumphi(17)=zsumphi(17)+znorm*HO%b_osc**2*(Input%ep*zQ21pt+Input%en*zQ21nt)
             zsumphi(18)=zsumphi(18)+znorm*HO%b_osc**2*(Input%ep*zQ22pt+Input%en*zQ22nt)
          else
             zsumphi(14)=zsumphi(14)+znorm*zQ2_2pt
             zsumphi(15)=zsumphi(15)+znorm*zQ2_1pt
             zsumphi(16)=zsumphi(16)+znorm*zQ20pt
             zsumphi(17)=zsumphi(17)+znorm*zQ21pt
             zsumphi(18)=zsumphi(18)+znorm*zQ22pt
          endif

        endif

!    ...............................
!    calculation of n-body density        
!    ...............................

!       print *,' computing Rho1B' 
        call ZRho1B_Angles(iang,ZRO,znorm,ZOBDJ0,NLEV) 
        if(Input%idens.ne.0) then
           call ZRho2B_Angles(iang,ZRO,ZKAPA10,ZKAPA01,znorm,NLEV)
           if(Input%iRho3B.eq.1) & 
           call ZRho3B_Angles(iang,ZRO,ZKAPA10,ZKAPA01,znorm,ZRO3BJJ,NLEV)
        endif
!     ................................................
        END DO !Fomenko
        END DO
!      ........................ end PNP
          call ZRho1B_Integral(iang,ZOBDJ0,0,zffJ0,ZRho1BJ0) 
        if(Input%idens.ne.0) then 
          call ZRho2B_Integral(iang,0,zffJ0)
          if(Input%iRho3B.eq.1) &
     &    call ZRho3B_Integral(iang,ZRO3BJJ,0,zffJ0,ZRho3BJJ)
        endif

        ! normalize all the quantities after PNP
        do ii=1,NOBS
           ZRES(iang,ii)=zsumphi(ii)/(PNP%NFOM**2)
        end do 

      ! output quantities after PNP but at different rotation Euler angle  
!      if(lpr) then
!       write(121,'(a,i4,a,3f8.3,a,17f12.5)')              &
!     &        '<phi1|R(',iang,':', alp, bet, gam, ')|phi2>=', &
!     &        dreal(ZRES(iang,1)),dreal(ZRES(iang,2)/ZRES(iang,1))+Nucl%ncore(0),   &
!     &        dreal(ZRES(iang,3)/ZRES(iang,1))+Nucl%ncore(1),dreal(ZRES(iang,7)/ZRES(iang,1)),&
!     &        dreal(ZRES(iang,8)/ZRES(iang,1)),ZRES(iang,10),ZRES(iang,11),ZRES(iang,12),ZRES(iang,15),ZRES(iang,16),ZRES(iang,17)      
!      endif

!       write(120,'(i5,3f8.3,7f12.5)')              &
!     &        iang, alp, bet, gam, &
!     &        dreal(ZRES(iang,1)),dreal(ZRES(iang,2)/ZRES(iang,1))+Nucl%ncore(0),   &
!     &        dreal(ZRES(iang,3)/ZRES(iang,1))+Nucl%ncore(1),dreal(ZRES(iang,7)/ZRES(iang,1)),&
!     &        dreal(ZRES(iang,8)/ZRES(iang,1)),ZRES(iang,16)

!      write(*,*) 'Q20p=',dreal(ZRES(iang,11)/ZRES(iang,1))
!     ................................................
        END DO  !end do loop beta
        END DO  !end do loop gamma
        END DO  !end do loop alpha



      ! < J1; q1,K1|| Q_lambda || J2; q2,K2> 
      ! external loop over angular momentum of initial state
!     .............................................
      !zqJ0(-2:2,2) = zzero

      Kernel%qpred(:,:,:,:) = zzero
      Kernel%qcpred(:,:,:,:) = zzero

!     .............................................
      do J2=0,JJmax,1
!     .............................................
          ZsumJKK(:,:,:) = zzero

!      ................ integration over angles
            ieuler=0
            !ztest = zzero
            DO ibet=1,max(1,NLEG_BET_Max) !start do loop in beta
               bet    = AMP%beta(ibet)
               zwb    = zone*AMP%wwbet(ibet)
               cosbet = cos(bet)

               if(AMP%NLEG_BET.eq.1) then
                  bet   = zero
                  zwb   = zone
                  cosbet = 1.d0
               endif

              DO ialp=1,AMP%NLEG_ALP !start do loop in alpha
                 alp=AMP%alpha(ialp)
                 zwa=AMP%wwalp(ialp)


              DO igam=1,AMP%NLEG_GAM !start do loop in gamma
                 gam=AMP%gamma(igam)
                 zwg=AMP%wwgam(igam)


                 ieuler=ieuler+1
                 if(ieuler.ne.ieuler_ang(ialp,ibet,igam)) &
                 stop 'ieuler_ang is not given properly'


                 iang1 = ieuler !ieuler_ang(ialp,ibet,igam)
                 ! weights of rotation angles: measure
                 zfac   = zwa*zwb*zwg  !*Ncosbet
                 zwdvol = zfac

                ! *******************************************  
                ! io: scalar (L=0) observables
                ! *******************************************  
                ! projection of J2(initial) onto z-direction
                do K2 =-J2,J2,1
                do K1 =-J2,J2,1
                 ! Wigner D-function
                   zap=cdexp(Zimag*K1*alp)
                   zbp=zone*dwignerI_gen('i',J2,K1,K2,bet)
                   zgp=cdexp(Zimag*K2*gam)
                   za1 = zap*zgp*zbp   
                   do io=1,8
                      zo1 = ZRES(iang1,io)   
                      ZsumJKK(io,K1,K2)=ZsumJKK(io,K1,K2) &
                        +weightJ(J2)*za1*zo1*zwdvol 
                   end do
                 !  if(K1.eq.K2 .and. K1.eq.0)  & 
                 !   ztest = ztest + weightJ(J2)*zap*zgp*zbp*zwdvol
                 end do ! K1
                 end do ! K2
                ! *******************************************  
                ! io: tensor (L=2, Delta_L=0) observables
                ! *******************************************  
                ! ZRES: <NZ,q1| O_2 |NZ,q2(alp,bet,gam)> 
                zqlm(-2,2) = ZRES(iang1,9) 
                zqlm(-1,2) = ZRES(iang1,10) 
                zqlm(0,2)  = ZRES(iang1,11) 
                zqlm(1,2)  = ZRES(iang1,12)
                zqlm(2,2)  = ZRES(iang1,13)

                zqlmc(-2,2) = dconjg(ZRES(iang1,14))
                zqlmc(-1,2) = dconjg(ZRES(iang1,15))
                zqlmc(0,2)  = dconjg(ZRES(iang1,16))
                zqlmc(1,2)  = dconjg(ZRES(iang1,17))
                zqlmc(2,2)  = dconjg(ZRES(iang1,18))

              ! if(J2.eq.2) then
              !   do lmu=-2,2  !mu 
              !      zqJ0(lmu,2) = zqJ0(lmu,2) + cdexp(Zimag*2*alp)*cdexp(Zimag*2*gam) &
              !       *dwignerI_gen('i',J2,2,2,bet)*zqlm(lmu,2)*zwdvol*weightJ(J2)
              !   enddo
              ! endif

              ! **************************
              ! spectroscopic Q_2 value (Jf=Ji)   <J1=J2,K1,q1 ||e Q2|| J2,K2,q2> 
              ! **************************
               J1 = J2 
               do K2 =-J2,J2,1
               do K1 =-J2,J2,1
               do lmu=-2,2  !mu 
                   K2p = K1-lmu !miis = -lmu
                   zap=cdexp(Zimag*K2p*alp)
                   zbp=zone*dwignerI_gen('i',J2,K2p,K2,bet)
                   zgp=cdexp(Zimag*K2*gam)

                   za1 = zap*zgp*zbp

                    Kernel%qpred(J2,K1,K2,0)            &
     &            = Kernel%qpred(J2,K1,K2,0)            &
     &              + weightJ(J2)*(2*J1+1.0)*iv(J1-K1)  &
     &              * wignei(J1,2,J2,-K1,lmu,K2p)       &
     &              * za1*zqlm(lmu,2)*zwdvol

               enddo !mlmu 
             enddo ! k1           
            enddo ! k2           

!     .........................................................   
!     io: tensor (L=2, Delta_L=2) observables
!     reduced matrix elements in BE2 value 
!     <J1 K1 q1|| Q2 || J2 K2 q2> =  <J+2,K1,q1 ||e Q2|| J,K2,q2> 
!     .........................................................   

            J1 = J2 + 2
            do K2 =-J2,J2,1
             do K1 =-J1,J1,1
               do mu=-2,2  !mu 
                   K2p = K1-mu !miis = -lmu
                   !if(abs(K2p).gt.J2) cycle


                   zap=cdexp(Zimag*K2p*alp)
                   zbp=zone*dwignerI_gen('i',J2,K2p,K2,bet)
                   zgp=cdexp(Zimag*K2*gam)

                   za1 = zap*zgp*zbp
                   Kernel%qpred(J2,K1,K2,2)            &
     &           = Kernel%qpred(J2,K1,K2,2)            &
     &              + weightJ(J2)*(2*J1+1.0)*iv(J1-K1) &
     &              * wignei(J1,2,J2,-K1,mu,K2p)       &
     &              *za1*zqlm(mu,2)*zwdvol

               enddo !mu 
             enddo ! k1           
            enddo ! k2

!     ......exchange q1 and q2 (Not K1 and K2): 

            do J1=J2,min(J2+2,JJmax)
              do K2 =-J2,J2,1
              do K1 =-J1,J1,1
              do mu=-2,2   
                   do nu=-2,2  !nu 
                      K2p = K1-mu !miis = -lmu
                      !if(abs(K2p).gt.J2) cycle

                   zam=cdexp(-Zimag*(K2+nu)*alp)
                   zgm=cdexp(-Zimag*(K2p+mu)*gam)
                   zbp=zone*dwignerI_gen('i',J2,K2,K2p,bet)
                   zbp2=zone*dwignerI_gen('i',2,nu,mu,bet)

                   za1 = zam*zgm*zbp*zbp2

                      Kernel%qcpred(J2,K1,K2,J1-J2)         &
     &             =  Kernel%qcpred(J2,K1,K2,J1-J2)         &
     &              + weightJ(J2)*(2*J1+1.0)*iv(abs(J1-K1)) &
     &              * wignei(J1,2,J2,-K1,mu,K2p)            &
     &              * za1*zqlmc(nu,2)*zwdvol
     !&              * za1*zqlm(nu,2)*zwdvol                 &

               enddo !nu 
               enddo !mu 
             enddo ! k1           
            enddo ! k2
            enddo ! J1 
         ! endif
!       end loops over Euler angles 
        END DO
        END DO !end do loop in beta
        END DO
        !print *,'jproj=',J2,ztest
!      ...........
!      ............................. print out density matrix elements for g.s.
        print *, '  '
        if(J2.eq.0) then
          call Write_ME1B(69,ZsumJKK(1,0,0),ZRho1BJ0) 
          !call Write_Rho1B4IMSRG(Rho1B_Full)     
          if(Input%idens.ne.0) then 
             call Write_ME2B(79,ZsumJKK(1,0,0),Rho1B_Full)
             if(Input%iRho3B.eq.1) & 
             call write_ME3B(89,ZsumJKK(1,0,0),ZRho3BJJ)
          endif
        endif


       do K1=-J2,J2
       do K2=-J2,J2 !ZsumJKK(io,K1,K2)
          !hjkk(J2,K1,K2) = zero 
          if(abs(ZsumJKK(1,K1,K2)).lt.CHOP) cycle
          Kernel%njkk(J2,K1,K2)    = ZsumJKK(1,K1,K2)
          Kernel%zz(J2,K1,K2)      = (ZsumJKK(2,K1,K2)/ZsumJKK(1,K1,K2)) + Nucl%ncore(1)
          Kernel%nn(J2,K1,K2)      = (ZsumJKK(3,K1,K2)/ZsumJKK(1,K1,K2)) + Nucl%ncore(0)
          Kernel%hjkk(J2,K1,K2)    = ZsumJKK(7,K1,K2)
          Kernel%J2(J2,K1,K2)      = ZsumJKK(8,K1,K2)/ZsumJKK(1,K1,K2)


         if(abs(ZsumJKK(1,K1,K2)).gt. 1.d-4) then 
         write(*,14) J2, K1,K2, dreal(ZsumJKK(1,K1,K2)),        &
     &            dreal(ZsumJKK(7,K1,K2)/ZsumJKK(1,K1,K2)),   &
     &            dreal(ZsumJKK(2,K1,K2)/ZsumJKK(1,K1,K2)),   &         
     &            dreal(ZsumJKK(3,K1,K2)/ZsumJKK(1,K1,K2)),   &
     &            dreal(ZsumJKK(8,K1,K2)/ZsumJKK(1,K1,K2)),   &
     &            Kernel%qpred(J2,K1,K2,2), &  !/ZsumJKK(1,K1,K2))
     &            Kernel%qcpred(J2,K1,K2,2)  !/ZsumJKK(1,K1,K2))
         endif
        enddo
        enddo
!     ----------- 
      end do ! end do loop in initial angular momentum J

       !    do lmu=-2,2  !mu 
       !       print *, lmu, zqJ0(lmu,2)
       !    enddo

       t1 = omp_get_wtime()
       print *, ' Time elapse:',t1-t0

      return


   14 format(' J=',i3,'(',2i3,') N^J=',f12.6,      &
     &       ' H^J=',f12.6,' Z(N)^J=',2f6.2,    &
     &       ' J(J+1)=',f8.3,' Qp=',2f12.8,' Qcp=',2f12.8)

       end
