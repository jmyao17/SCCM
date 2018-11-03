        subroutine Kernel4Axial(lpr,Nphi1,Nphi2,NLEV)
!       ......................................
!       If mphi=1: no PNP
!       If pairing collapse: no PNP
!       ......................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)
        logical lpr
        Integer a12,a34,b12,Ncosbet                                     ! symmetry for beta angle 
        Integer ieuler_ang(1:AMP%NLEG_ALP,1:AMP%NLEG_BET,1:AMP%NLEG_GAM)
        real*8 start1,finish1,start,finish
        REAL(DP) CG

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
        !DIMENSION ZLambda2B (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        Complex*16, Allocatable, DIMENSION(:,:) :: ZRO3BJJ,ZRho3BJJ,ZLambda3B    

        REAL*8, DIMENSION(0:1,1:HO%nMax,0:HO%lmax*2,0:1,1:HO%nMax,0:HO%LMax*2) :: Rho1B_Full 

        complex*16, dimension(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2) :: ZTD,ZTDc
        complex*16, dimension(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2) :: ZTDF,ZTDFc  

        !MATRIX ELEMENTS
        DIMENSION ZROT(NLEV,NLEV)
        DIMENSION ZAM_phi_omega(NLEV,NLEV)

        !OTHER QUANTITIES

        DIMENSION ZZc0(NLEV,NLEV),ZZc1(NLEV,NLEV)
        DIMENSION Zsumphi(NOBS)
        DIMENSION ZRES(AMP%NNNP,NOBS)
        DIMENSION ZsumJJ(NOBS)
        CHARACTER*100 cwf,cwf1,cwf2
        complex*16 zqlm(-2:2,2),zqlm2(-2:2,2)

        real*8 t0,t1,t2,t3,t4,t5,t6,t7,t8 !,omp_get_wtime()

        write(*,*) ' -> Start to compute kernels ...'
!  ....................... preparation
        Aneut = Input%nneut !Nucl%nucleon(0) - Nucl%ncore(0)
        Aprot = Input%nprot !Nucl%nucleon(1)
        aaa=1.d0*NLEV*(NLEV+1)/2.d0
        snpfaf=(-1.d0)**(aaa)         !for the Pfaffian, check it!
        NLEG_BET_Max = AMP%NLEG_BET

!  ....................... Use symmetry (2), otherwise, (1)
        Ncosbet = 2         ! beta ranges from 0 to pi/2
        if(Ncosbet.eq.2 .and. AMP%NLEG_BET/2.le.0) NLEG_BET_Max=2
!  ..........................
        if(Ncosbet.eq.0) stop 'Ncosbet should be 1 or 2'
        do J=0,JJmax
           weightJ(J) = (2.d0*J+1.d0)/2.d0  !(8.d0*pi**2)  
        enddo
!  read wave functions
!  ................................ read u0,v0 of ref. state 1
      if(.NOT. ALLOCATED(HFB%U0))  ALLOCATE(HFB%U0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%V0))  ALLOCATE(HFB%V0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%U1))  ALLOCATE(HFB%U1(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%V1))  ALLOCATE(HFB%V1(1:HO%NLEV,1:HO%NLEV))
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


!    ................................ read u1,v1 of ref. state 2
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
!     ..................... check if anyone of the two wavefunctions is HF state
!     ide(0:1) = 0 (HF); 1 (HFB)
      do it=0,1
         HFB%ide(it)=HFB%ide0(it)*HFB%ide1(it)
      enddo !it


      if(.NOT. ALLOCATED(AMP%ZROT_m1m2))  ALLOCATE(AMP%ZROT_m1m2(1:HO%NLEV,1:HO%NLEV))
!       ................ determine u and v in canonic basis
!         call canon(NLEV)
!       ................ truncation canonic basis 
!         call dimens(NLEV)
       
          print *,' --> writing one-body transition density into:',File%TD1B
            open(59,file=File%TD1B,status='unknown')
          print *,' --> writing one-body density into:',File%Rho1B
            open(69,file=File%Rho1B,status='unknown')
         if(Input%idens.ne.0) then 
          print *,' --> writing two-body density into:',File%Rho2B
            open(300,file=File%Rho2B,status='unknown')
            if(Input%iRho3B.eq.1) &
            open(89,file=File%Rho3B,status='unknown')
         endif


!        ......... for shell-model space, alwasy using PNP
         if(Input%IntType .eq. 0) then
            HFB%ide(0) = 1 
            HFB%ide(1) = 1
         endif
!        if there is pairing in both neutrons and protons, use
!        pfaffian method 
! 
!        determine the normalization factor
         if(HFB%ide(0)+HFB%ide(1) .ge.1) then  
            ! using Neergard-Wuest method
            call zbareover(ZV0,ZU0,ZZc0,Det0,NLEV)  !q (bra)
!        ......... 10^N  * 10^N is from pfaffian
            zDetvac0=Det0*10.d0**(-NLEV*Input%iscale/2) !a bug is removed here, zone*((1./Det0)**.25)
            call zbareover(ZV1,ZU1,ZZc1,Det1,NLEV)  !q' (ket)
            zDetvac1=Det1*10.d0**(-NLEV*Input%iscale/2) !zone*((1./Det1)**.25)
            write(*,*) '<phi|phi>(NW):',zDetvac0,zDetvac1
            ! using pfaffian method
            !call normalization(ZV0,ZU0,ZNF,NLEV) 
            !write(*,*) '<phi_0|phi_0>(pf):',ZNF
            !call Pfaf_Over(ZZc0,ZU0,ZV0,ZdetPfaf,NLEV)
            !write(*,*) '<phi|phi>(pf):',ZdetPfaf
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
      ZTDF(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)=zzero
      ZTDFc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)=zzero

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
     
            t0 = omp_get_wtime()
            Dens%ZRho2BJJ(:,:,:) = zzero
            t1 = omp_get_wtime()
            print *, 'initializing Dens%ZRho2BJJ:',t1-t0
       endif
!  .......................................
        iang=0
! .......................... loop over beta angle
!        symmetry
        DO k=1,max(NLEG_BET_Max/Ncosbet,1)
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
           ieuler_ang(kkk,k,kk) = iang
           iang_pimbet = iang+(AMP%NLEG_BET-2*k+1)
           zffJ0=zwa*zwb*zwg*Ncosbet ! only for the DMEs of g.s. (J=0) 

!       ........... TWO-BODY Density coupled to J=0
!        The two-body density requires large memory, may cause
!        fragmentation error

            !print *, ' Initializing one-body density ...'
            ZOBDJ0 (0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)       = zzero
            ZTD(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)   = zzero 
            ZTDc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)   = zzero 


            !print *, ' ->>FSPB%tnlj_fmax=',FSPB%tnlj_fmax
        if(Input%idens.ne.0) then
            t0 = omp_get_wtime()
            !print *, ' Initializing two-body density ...'
            Dens%ZTBDJJ (:,:,:) = zzero          
            if(Input%iRho3B.eq.1)  then 
            !print *, ' Initializing three-body density ...'
            ZRO3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax)  = zzero
            endif
            t1 = omp_get_wtime()
            print *, ' Time consummed for initializing densities:',t1-t0
        endif
!+-------------------------------------------------------------------------c
!       Rotation matrix                                                     c
!+-------------------------------------------------------------------------c
       call ROTATION_MAT(alp,bet,gam,ZROT,NLEV)
       do ia=1,NLEV
       do ib=1,NLEV
          AMP%ZROT_m1m2(ia,ib) = ZROT(ia,ib)
       enddo
       enddo
!+-------------------------------------------------------------------------c
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

          ! write(*,'(a,i4,i4,a)') '  Gauge angles: (',LLLP,LLLN,')'
!          call cpu_time(start)

          phiLN=(2*pi/PNP%NFOM)*LLLN               ! phi_n
          if(PNP%NFOM.eq.1) phiLN=0.d0

          zphas_n=cdexp(-Zimag*phiLN*Aneut)     ! exp(-i*N*phi_n)

          PNP%phi(1) = phiLN
          PNP%phi(2) = phiLP
          PNP%zemiphi(1) = zphas_n
          PNP%zemiphi(2) = zphas_p
!         ...............
!         rotated U and V
!         ...............
!         print *, zphas_n,zphas_p
!        if(LLLN.eq.0 .and. LLLP.eq.0) t1 = omp_get_wtime()

         call UVphiomega(phiLP,phiLN,ZROT,ZU1,ZV1,ZU1bar,ZV1bar,ZAM_phi_omega,NLEV)

!        if(LLLN.eq.0 .and. LLLP.eq.0) t2 = omp_get_wtime()
!        if(LLLN.eq.0 .and. LLLP.eq.0) print *, 'UVphi:',t2-t1 

        !norm overlap   
!        ........ Onishi formula works
         if(PNP%NFOM .eq.1) then
         call Onishi(ZU0,ZV0,ZU1bar,ZV1bar,znorm,HO%NLEV)
!          write(*,*) ' Method I Onishi formula: znorm=',znorm         

!        ........ Neegard_Wust subroutine needs further check
!         call Neegard_Wust(ZU0,ZV0,ZU1bar,ZV1bar,ZZTZ,NLEV,znorm)
!         call zdiag(NLEV,ZZTZ,znorm)
!         znorm = znorm*zDetvac0*zDetvac1*zphas_p*zphas_n
!         write(*,*) ' Method I: znorm=',znorm !,zDetvac0,zDetvac1
         else
           call Pfaf_Over(ZZc0,ZU1bar,ZV1bar,ZdetPfaf,NLEV)
           znorm = snpfaf*zDetvac0*zDetvac1*ZdetPfaf*zphas_p*zphas_n
!        write(*,*) ' Method II: znorm=',znorm !,zDetvac0*zDetvac1,ZdetPfaf
         endif

        zsumphi(1)=zsumphi(1)+znorm

        call ROKAPPAphiomega(ZAM_phi_omega,ZU0,ZV0,ZU1,ZV1,ZRO,ZKAPA10,ZKAPA01,NLEV)
        call ZN_0(ZRO,Zprot_0,Zneut_0,NLEV)

        zsumphi(2)=zsumphi(2)+znorm*Zprot_0
        zsumphi(3)=zsumphi(3)+znorm*Zneut_0

!  .......... two ways to calculate the energy overlap
!  (1)
!       call cpu_time(start)
!        call system_clock ( nt3, nclock_rate, nclock_max )
!        call HFB_FIELD_COMPLEX(ZRO,zkapa10,zkapa01,H%zME1BM,zgamma,zham,zdelta10,zdelta01,NLEV)  
!        call HFB_ENER_comp(zro,zkapa01,H%zME1BM,zgamma,zdelta10,zEkin,zEHF,zEPa,zEHFB,NLEV)
!        write(*,*) 'Method I: EHFB=',zEkin,zEHFB
!        call system_clock ( nt4, nclock_rate, nclock_max )
!        write ( *, * ) 'Elapsed real time = ', real ( nt4 - nt3 ) / real ( nclock_rate )
!       call cpu_time(finish)
!       print '("Time = ",f10.3," seconds.")',finish-start
!  (2)

          call HFB_Energy_COMPLEX(zro,zkapa10,zkapa01,H%zME1BM,zEkin,zEHFB,NLEV,.false.)

!        if(LLLN.eq.0 .and. LLLP.eq.0) t5 = omp_get_wtime!()
!        if(LLLN.eq.0 .and. LLLP.eq.0) print *, 'Energy_overlap:',t5-t4 
!         write(*,*) 'Method II: EHFB=',zEkin,zEHFB
!  ........................................................................
        zsumphi(4)=zsumphi(4)+znorm*zEkin
        zsumphi(5)=zsumphi(5)+znorm*zEHF
        zsumphi(6)=zsumphi(6)+znorm*zEPa
        zsumphi(7)=zsumphi(7)+znorm*(zEHFB+H%E0)
!        write(*,*) ' Start ZJ2_0 ....'
        call ZJ2_0(ZRO,ZKAPA10,ZKAPA01,zJ2_MV,NLEV)

!        if(LLLN.eq.0 .and. LLLP.eq.0) t6 = omp_get_wtime()
!        if(LLLN.eq.0 .and. LLLP.eq.0) print *, '<J2>:',t6-t5 

        zsumphi(8)=zsumphi(8)+znorm*zJ2_MV
!   ............................... quadrupole moments

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
           zsumphi(11)=zsumphi(11)+znorm*HO%b_osc**2*(Input%ep*zQ20p +Input%en*zQ20n)
           zsumphi(12)=zsumphi(12)+znorm*HO%b_osc**2*(Input%ep*zQ21p +Input%en*zQ21n)
           zsumphi(13)=zsumphi(13)+znorm*HO%b_osc**2*(Input%ep*zQ22p +Input%en*zQ22n)

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
             zsumphi(15)=zsumphi(15)+znorm*HO%b_osc**2*(Input%ep*zQ2_1pt+Input%en*zQ2_1n)
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
!      ...............................
          call ZQ1B(cME1B%r2*zone,ZRO,zQ0p,zQ0n,NLEV)
          if(Input%IntType.eq.0) then
            zsumphi(19)=zsumphi(19)+znorm*(Input%ep*zQ0p+Input%en*ZQ0n)
          else
            zsumphi(19)=zsumphi(19)+znorm*zQ0p
          endif
!    ...............................
!    calculation of n-body density        
!    ...............................
        call ZRho1B_Angles(iang,ZRO,znorm,ZOBDJ0,NLEV) 
        !print *,' computing TD1B' 
        call ZTD1B_Angles(iang,ZRO,znorm,ZTD,ZTDc,NLEV) 
        if(Input%idens.ne.0) then
!        print *,' computing Rho2B' 
           call ZRho2B_Angles(iang,ZRO,ZKAPA10,ZKAPA01,znorm,NLEV)
           if(Input%iRho3B.eq.1) & 
           call ZRho3B_Angles(iang,ZRO,ZKAPA10,ZKAPA01,znorm,ZRO3BJJ,NLEV)
        endif
!     ................................................
        END DO !Fomenko
        END DO
!      ........................ end PNP
        call ZRho1B_Integral(iang,ZOBDJ0,0,zffJ0,ZRho1BJ0) 
!       axial case, Ki=Kf=0
        call ZTD1B_Integral(alp,bet,gam,ZTD,ZTDc,zffJ0,ZTDF,ZTDFc) 
        if(Input%idens.ne.0) then 
          call ZRho2B_Integral(iang,0,zffJ0)
          if(Input%iRho3B.eq.1) &
     &    call ZRho3B_Integral(iang,ZRO3BJJ,0,zffJ0,ZRho3BJJ)
        endif

        do ii=1,NOBS
           ZRES(iang,ii)=zsumphi(ii)/(PNP%NFOM**2)
        end do
!   ......... timing

!      print *, 'UVphiomega:',(t2-t1)*PNP%NFOM**2
!      print *, 'Pfaffian:',(t3-t2)*PNP%NFOM**2
!      print *, 'ROKAPPAphiomega:',(t4-t3)*PNP%NFOM**2
!      print *, 'EHFB:',(t5-t4)*PNP%NFOM**2
!      print *, 'zJ0:',(t6-t5)*PNP%NFOM**2
!      print *, 'Calc.Dens:',(t7-t6)*PNP%NFOM**2
!      print *, 'Sum_of_Time:',t8-t0

       write(*,'(a,i3,a,3f8.3,a,7f12.5)')              &
     &        '<phi1|R(',iang,':', alp, bet, gam, ')|phi2>=', &
     &        dreal(ZRES(iang,1)),dreal(ZRES(iang,2)/ZRES(iang,1))+Nucl%ncore(0),   &
     &        dreal(ZRES(iang,3)/ZRES(iang,1))+Nucl%ncore(1),dreal(ZRES(iang,7)/ZRES(iang,1)),&
     &        dreal(ZRES(iang,8)/ZRES(iang,1)),ZRES(iang,16)     

!       write(120,'(i5,3f8.3,5f12.5)')              &
!     &        iang, alp, bet, gam, &
!     &        dreal(ZRES(iang,1)),dreal(ZRES(iang,2)/ZRES(iang,1))+Nucl%ncore(0),   &
!     &        dreal(ZRES(iang,3)/ZRES(iang,1))+Nucl%ncore(1),dreal(ZRES(iang,7)/ZRES(iang,1)),&
!     &        dreal(ZRES(iang,8)/ZRES(iang,1))
!      write(*,*) 'Q20p=',dreal(ZRES(iang,11)/ZRES(iang,1))
!     ................................................
        END DO  !end do loop beta
        END DO  !end do loop gamma
        END DO  !end do loop alpha
!      ........................ end AMP 


!     start do loop in initial angular momentum Ji

      DO J=0,JJmax,2

         iis = J
         K   = 0
         Kp  = 0
!......... initialization
         Kernel%qpred(iis,0,0,0:2) = zzero 
         Kernel%qcpred(iis,0,0,0:2) = zzero 
         Kernel%q0p(iis,0,0) = zzero 
         Kernel%q0pc(iis,0,0) = zzero 


!      .........................................
            do kk=1,NOBS
               ZsumJJ(kk)=Zzero
            end do

!      ................ integration over angles

        ieuler=0
        test = zero
        DO ibet=1,max(1,NLEG_BET_Max/Ncosbet) !start do loop in beta

         bet=AMP%beta(ibet)
         zwb=zone*AMP%wwbet(ibet)
         cosbet = cos(bet)
         if(AMP%NLEG_BET.eq.1) then
            bet   = zero
            zwb   = zone
            cosbet = 1.d0 
         endif

         zbp=zone*dwignerI_gen('i',J,K,Kp,bet)
         zbm=zone*dwignerI_gen('i',J,K,-Kp,bet)

        DO ialp=1,AMP%NLEG_ALP !start do loop in alpha
         alp=AMP%alpha(ialp)
         zwa=AMP%wwalp(ialp)

         zap=cdexp(Zimag*K*alp)
         zam=cdexp(-Zimag*K*alp)

        DO igam=1,AMP%NLEG_GAM !start do loop in gamma
         gam=AMP%gamma(igam)
         zwg=AMP%wwgam(igam)

         zgp=cdexp(Zimag*Kp*gam)
         zgm=cdexp(-Zimag*Kp*gam)

!......... symmetry
           ieuler=ieuler+1

         if(ieuler.ne.ieuler_ang(ialp,ibet,igam)) &
     &   stop 'ieuler_ang is not given properly'


         zff  = zwa*zwb*zwg*(zap*zbp*zgp)*Ncosbet
         zfac = zwa*zwb*zwg*Ncosbet

!  .....................................
         do jj=1,8 !NOBS
            ZsumJJ(jj)=ZsumJJ(jj)+weightJ(J)*zff*ZRES(ieuler,jj)
         end do
         zqlm(-2,2) = ZRES(ieuler,9) 
         zqlm(-1,2) = ZRES(ieuler,10) 
         zqlm( 0,2) = ZRES(ieuler,11) 
         zqlm( 1,2) = ZRES(ieuler,12) 
         zqlm( 2,2) = ZRES(ieuler,13) 
         zqlm2(-2,2) = ZRES(ieuler,14)
         zqlm2(-1,2) = ZRES(ieuler,15)
         zqlm2( 0,2) = ZRES(ieuler,16)
         zqlm2( 1,2) = ZRES(ieuler,17)
         zqlm2( 2,2) = ZRES(ieuler,18)

         zq00p = ZRES(ieuler,19)


!     .........................................................  <J,q1||er^2|| J,q2>
           Kernel%q0p(iis,0,0) = Kernel%q0p(iis,0,0) &
     &              + weightJ(iis)*zq00p             &
     &              * dwignerI_gen('i',iis,0,0,bet)  & !sdjmk(iis,iis,iis,ibet,0)
     &              * zfac

!     .........................................................  <J,q2||er^2|| J,q1>
           Kernel%q0pc(iis,0,0) = Kernel%q0pc(iis,0,0) &
     &              + weightJ(iis)*DCONJG(zq00p)       &
     &              * dwignerI_gen('i',iis,0,0,bet)    & !sdjmk(iis,iis,iis,ibet,0)
     &              * zfac 

!     ......................................................... spectroscopic Q_2 value 
            ifs = iis
            do lmu=-2,2  !mu 

!     .........................................................  <J,q1 ||e Q2|| J,q2> 
               miis = -lmu

               Kernel%qpred(iis,0,0,0) = Kernel%qpred(iis,0,0,0)      &
     &              + weightJ(iis)*(2*ifs+1)*iv(ifs)*zqlm(lmu,2)& ! *on
     &              * dwignerI_gen('i',iis,miis,0,bet)            &
     &              * wignei(iis,2,ifs,miis,lmu,0)                &
     &              * zfac 
!     .........................................................  <J,q2 ||e Q2|| J,q1>
             ! if(Nphi1.ne.Nphi2) then
              do nu = -2, 2
               Kernel%qcpred(iis,0,0,0) = Kernel%qcpred(iis,0,0,0)      &
     &              + weightJ(iis)*(2*ifs+1)*iv(ifs)*iv(abs(lmu-nu)) &
     &              * dwignerI_gen('i',2,-nu,-lmu,bet)                 &
     &              * DCONJG(zqlm2(nu,2))                               &
     &              * dwignerI_gen('i',iis,0,miis,bet)                 &
     &              * wignei(iis,2,ifs,miis,lmu,0)                     &
     &              * zfac 
              enddo ! nu
             ! endif
           enddo !mlmu 
           
!     ......................................................... reduced matrix elements in BE2 value 
            ifs = iis + 2
            do lmu=-2,2  !mu   
!     .........................................................   
!     <Jf qf|| Q2 || Ji qi> =  <J+2,q1 ||e Q2|| J,q2> 
!     .........................................................   
               miis = -lmu
               Kernel%qpred(iis,0,0,2) = Kernel%qpred(iis,0,0,2)      &
     &              + weightJ(iis)*sqrt(2*ifs+1.0)*zqlm(lmu,2)  &  ! *on
     &              * CG(iis*2,2*miis,4,2*lmu,ifs*2,0)          &
     &              * dwignerI_gen('i',iis,miis,0,bet)       &
     &              * zfac 
!   ...........way two
!     &              + weightJ(iis)*(2*ifs+1)*iv(ifs)*zqlm(lmu,2)  &  ! *on
!     &              * dwignerI_gen('i',iis,miis,0,bet)            & ! djmk(iis,miis,0,cosbet,0)*iv(miis) 
!     &              * wignei(iis,2,ifs,miis,lmu,0)                &  
!     &              * zfac !*sqrt(16*pi/5.d0)
!     ......................................  exchange q1 and q2 with symmetry
              !if(Nphi1.ne.Nphi2) then
              do nu = -2, 2
               Kernel%qcpred(iis,0,0,2) = Kernel%qcpred(iis,0,0,2)         &
     &              + weightJ(iis)*(2*ifs+1)*iv(ifs)*iv(abs(lmu-nu))   &
     &              * dwignerI_gen('i',2,-nu,-lmu,bet)                 &
     &              * DCONJG(zqlm2(nu,2))                              &
     &              * dwignerI_gen('i',iis,0,miis,bet)                 &
     &              * wignei(iis,2,ifs,miis,lmu,0)                     &
     &              * zfac !*sqrt(16*pi/5.d0)
                enddo ! nu
             ! endif
            enddo !mlmu 
        END DO
        END DO !end do loop in beta
        END DO

!      ............................. print out density matrix elements for g.s.

        if(J.eq.0) then
          call Write_TD1B(59,ZTDF,ZTDFc,zsumJJ(1)) 
          call Write_ME1B(69,zsumJJ(1),ZRho1BJ0) 
          !call Write_Rho1B4IMSRG(Rho1B_Full)     
          if(Input%idens.ne.0) then 
             call Write_ME2B(79,ZsumJJ(1),Rho1B_Full)
             if(Input%iRho3B.eq.1) & 
             call write_ME3B(89,ZsumJJ(1),ZRho3BJJ)
          endif
        endif



        OverNo=  dreal(ZsumJJ(1))
        OverProt=dreal(ZsumJJ(2))
        OverNeut=dreal(ZsumJJ(3))
        OverHFB=dreal(ZsumJJ(7))
        OverJ2=dreal(ZsumJJ(8))


!        END DO !end do loop Kprime
!        END DO !end do loop K
!      ..............................................
      Kernel%nn(iis,0,0)      = (ZsumJJ(3)/ZsumJJ(1)) + Nucl%ncore(0)
      Kernel%zz(iis,0,0)      = (ZsumJJ(2)/ZsumJJ(1)) + Nucl%ncore(1)
      Kernel%hjkk(iis,0,0)    = (ZsumJJ(7))
      Kernel%njkk(iis,0,0)    = (ZsumJJ(1))
      Kernel%J2(iis,0,0)      = (ZsumJJ(8)/ZsumJJ(1))
       
      write(*,14) iis, dreal(ZsumJJ(1)),        &
     &            dreal(ZsumJJ(7)/ZsumJJ(1)),   &
     &            dreal(ZsumJJ(2)/ZsumJJ(1)),   &         
     &            dreal(ZsumJJ(3)/ZsumJJ(1)),   &
     &            dreal(ZsumJJ(8)/ZsumJJ(1)),   &
     &            Kernel%qpred(iis,0,0,2), & !/ZsumJJ(1))
     &            Kernel%qpred(iis,0,0,2) !/ZsumJJ(1))

       end do !end do loop in initial angular momentum J
   14 format(' J=',i4,' N^J_00=',f10.4,              &
     &       ' H^J_00=',f10.4,' Z(N)^J_00=',2f6.2,   &
     &       ' J(J+1)=',f8.3,' Qp=',2f12.8,' Qpc=',2f12.8)

       end
