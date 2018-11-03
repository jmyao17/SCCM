!      Subroutines for computing two-body densities
!      .............................................
!      If found error while calling this subroutine
!      probably, the tnljm%level is not defined in Hamiltonian.f90       
!      .............................................
        subroutine ZRho2B_Angles(iangle,ZRO,ZKAPA10,ZKAPA01,znorm,NLEV)
        USE VAPHFB_PAR
        implicit none
        integer   iangle,NLEV 
        integer   k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 ZRO,ZKAPA10,ZKAPA01,Znorm,Znorm2
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION ZKAPA10(NLEV,NLEV)    ! density matrix
        DIMENSION ZKAPA01(NLEV,NLEV)    ! density matrix

        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER a_ij,b_ij,a_kl,b_kl
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l

        REAL(DP) CG
        real*8 omp_get_wtime,t0,t1,cg1,cg2
        complex*16 Trace
        real*8  Metric
!     ............. print out the one-body density
        t0 = omp_get_wtime()
        write(*,'(a30)',advance='no') '  Computing the Rho2B ....'
        Znorm2 = Znorm/(PNP%NFOM**2)

        Trace=zzero
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Trace,Znorm2,zkapa01,zro,zkapa10,SPB,TPB,tnljm,CG_Save,Dens) 
!$OMP DO schedule(dynamic)  REDUCTION(+:Trace) 
       do bm=0, TPB%bMax-1
           tt= TPB%btt(bm)
           PP= TPB%bPP(bm)
           JJ= TPB%bJJ(bm)
           bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        do b_ij=0, TPB%aMax(bb)-1
!-------------------------------------------------------
           a_ij= b_ij +1
           if(a_ij.gt.aMaxMax) stop 'aMaxMax is too small !'
           t_i = TPB%t1(bb,a_ij)
           n_i = TPB%n1(bb,a_ij)  ! 
           l_i = TPB%twol1(bb,a_ij)  ! l_i is doubled
           j_i = TPB%j1(bb,a_ij)  ! j_i is doubled
           LJ_i = (l_i+j_i-1)/2

           t_j = TPB%t2(bb,a_ij)
           n_j = TPB%n2(bb,a_ij)  ! 
           l_j = TPB%twol2(bb,a_ij)  ! l_j is doubled
           j_j = TPB%j2(bb,a_ij)  ! j_j is doubled
           LJ_j =(l_j+j_j-1)/2
        do b_kl=0, TPB%aMax(bb)-1
           a_kl= b_kl  +1
           if(a_kl.gt.aMaxMax) stop 'aMaxMax is too small !'
           t_k = TPB%t1(bb,a_kl)
           n_k = TPB%n1(bb,a_kl)  ! 
           l_k = TPB%twol1(bb,a_kl)  ! l_i is doubled
           j_k = TPB%j1(bb,a_kl)  ! j_i is doubled
           LJ_k = (l_k+j_k-1)/2

           t_l = TPB%t2(bb,a_kl)
           n_l = TPB%n2(bb,a_kl)  ! n2 starts from 1 
           l_l = TPB%twol2(bb,a_kl)  ! l_j is doubled
           j_l = TPB%j2(bb,a_kl)  ! j_j is doubled
           LJ_l = (l_l+j_l-1)/2

           if(t_i.eq.t_j.and.n_i.eq.n_j.and.LJ_i.eq.LJ_j) then
             Metric=(1+(-1)**(JJ))/2
           else
             Metric=2
           endif
!      .......... the n in SPB%n starts from 0, while n in others starts from 1 !!!
       if(SPB%VType(lj_i,n_i-1,t_i) .ne. 1) cycle
       if(SPB%VType(lj_j,n_j-1,t_j) .ne. 1) cycle
       if(SPB%VType(lj_k,n_k-1,t_k) .ne. 1) cycle
       if(SPB%VType(lj_l,n_l-1,t_l) .ne. 1) cycle

!      a bug was removed here: tnljm%level was not defined
!      .............. loops over m

       do m_i=-j_i,j_i,2
          k_i = tnljm%level(m_i,lj_i,n_i,t_i)
       do m_j=-j_j,j_j,2

          MM  = (m_i+m_j)/2

          k_j = tnljm%level(m_j,lj_j,n_j,t_j)
!          cg1 = CG(j_i,m_i,j_j,m_j,2*JJ,2*MM)
          cg1 = CG_Save((j_i+1)/2,(m_i+1)/2,(j_j+1)/2,(m_j+1)/2,JJ,MM) 
       do m_k=-j_k,j_k,2
          m_l = 2*MM - m_k

!         ................... a bug is removed here
          if(abs(m_l).gt.j_l) cycle

          k_k = tnljm%level(m_k,lj_k,n_k,t_k)
          k_l = tnljm%level(m_l,lj_l,n_l,t_l)
!          cg2=CG(j_k,m_k,j_l,m_l,2*JJ,2*MM)
          cg2 = CG_Save((j_k+1)/2,(m_k+1)/2,(j_l+1)/2,(m_l+1)/2,JJ,MM) 

!    .......................
!      <c^+_i c^+_j c_l c_k> 
!    = <c^+_i c_k  ><c^+_j c_l>
!    - <c^+_i c_l  ><c^+_j c_k>
!    - <c^+_i c^+_j><c_l c_k>
!    .......................

          Dens%ZTBDJJ(a_ij,a_kl,bb) = Dens%ZTBDJJ(a_ij,a_kl,bb)      &
     &         + cg1*cg2*Znorm2/(2.d0*JJ+1.d0)                        &
     &         *(ZRO(k_k,k_i)*ZRO(k_l,k_j)                   &
     &          -ZRO(k_l,k_i)*ZRO(k_k,k_j)                   & 
     &          +ZKAPA01(k_i,k_j)*ZKAPA10(k_k,k_l)) 
         enddo
         enddo
         enddo
        
           if(a_ij.eq.a_kl) Trace = Trace + (2.d0*JJ+1.d0)*Dens%ZTBDJJ(a_ij,a_kl,bb)*Metric

        enddo
        enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
        write(*,'(a5)',advance='no') 'done'
        t1 = omp_get_wtime()
        write(*,'(a20,f10.5)',advance='no'), '...Time elapsed:',t1-t0
        write(*,'(a14,2f15.8)') ' Trace[Rho2B]=',Trace !/Znorm2
        end


        subroutine ZRho2B_Integral(iangle,J,zff)
        USE VAPHFB_PAR
        implicit none
        integer   iangle,J
        integer   bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        complex*16 zff
        do bm=0, TPB%bMax-1
           tt= TPB%btt(bm)
           PP= TPB%bPP(bm)
           JJ= TPB%bJJ(bm)
           bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        do b_ij=0, TPB%aMax(bb)-1
           a_ij= b_ij +1
        do b_kl=0, TPB%aMax(bb)-1
           a_kl= b_kl  +1
           Dens%ZRho2BJJ(a_ij,a_kl,bb) = Dens%ZRho2BJJ(a_ij,a_kl,bb)  &
     &                           + weightJ(J)*zff*Dens%ZTBDJJ(a_ij,a_kl,bb)  

         enddo
         enddo
         enddo
       end

        subroutine Write_ME2B(ie,znorm,OBD_Save)
        USE VAPHFB_PAR
        implicit none
        integer   ie
        integer   J,k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        integer   bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        complex*16 zff,znorm,Rho2BN,Trace(0:1)
        REAL*8, DIMENSION(0:1,1:HO%nMax,0:HO%LMax*2,0:1,1:HO%nMax,0:HO%LMax*2) :: OBD_Save
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
        real*8  Metric
        integer  d_12,d_34,iphase
        integer  d_ik,d_jl,d_il,d_jk
        real*8  rho_ik,rho_jl,rho_il,rho_jk
        complex*16 dme,TBD
!     ...............................

 20   format(/)
 50   format('*** 2b DME from GCM calculation ***')
 60   format('# basis (block, aMax)')
195   format('bMax=',i6)
196   format(2i10)
197   format('# block',i6,':',' aMax=',i6)

      open(ie,file='Rho2B.dat',status='unknown')
      write(300,*)' a_ij   a_kl   bb   Rho2B'
      write(300,*)' -----------------------------------'
      write(ie,50)
      write(ie,60)
      write(ie,195) TPB%bMax-1
      if(TPB%bMax.gt.bMax) stop 'bMax is too small'
      do bm=0, TPB%bMax-1          ! The block starts from 0 to bMax_calc-1 
         write(ie,196) bm,TPB%cMax(bm) ! right 
      enddo

      write(ie,20)

       Trace(0:1) = zzero
       do bm=0, TPB%bMax-1
           write(ie,197) bm,TPB%cMax(bm)
           tt= TPB%btt(bm)
           PP= TPB%bPP(bm)
           JJ= TPB%bJJ(bm)
           bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        do b_ij=0, TPB%aMax(bb)-1
           a_ij= b_ij +1
        t_i = TPB%t1(bb,a_ij)
        n_i = TPB%n1(bb,a_ij)  ! 
        l_i = TPB%twol1(bb,a_ij)  ! l_i is doubled
        j_i = TPB%j1(bb,a_ij)  ! j_i is doubled
        LJ_i = (l_i+j_i-1)/2

        t_j = TPB%t2(bb,a_ij)
        n_j = TPB%n2(bb,a_ij)  ! 
        l_j = TPB%twol2(bb,a_ij)  ! l_j is doubled
        j_j = TPB%j2(bb,a_ij)  ! j_j is doubled
        LJ_j =(l_j+j_j-1)/2

       if(a_ij .ne. TPB%a(bb,FSPB%tnlj(lj_i,n_i,t_i),FSPB%tnlj(lj_j,n_j,t_j))) &
     &     stop 'Error in QNs for TPB%a'

        do b_kl=0, TPB%aMax(bb)-1
           a_kl= b_kl  +1
        t_k = TPB%t1(bb,a_kl)
        n_k = TPB%n1(bb,a_kl)  ! 
        l_k = TPB%twol1(bb,a_kl)  ! l_i is doubled
        j_k = TPB%j1(bb,a_kl)  ! j_i is doubled
        LJ_k = (l_k+j_k-1)/2

        t_l = TPB%t2(bb,a_kl)
        n_l = TPB%n2(bb,a_kl)  ! 
        l_l = TPB%twol2(bb,a_kl)  ! l_j is doubled
        j_l = TPB%j2(bb,a_kl)  ! j_j is doubled
        LJ_l = (l_l+j_l-1)/2

       if(t_i.eq.t_j.and.n_i.eq.n_j.and.LJ_i.eq.LJ_j) then
             Metric=(1+(-1)**(JJ))/2
        else
             Metric=2
        endif

!     ............................................................
        d_12 = 0
        d_34 = 0
        if(n_i.eq.n_j.and.LJ_i.eq.LJ_j) d_12=1
        if(n_k.eq.n_l.and.LJ_k.eq.LJ_l) d_34=1

!      factor = dsqrt((one+d_12)*(one+d_34))

        rho_ik = OBD_Save(t_i,n_i,lj_i,t_k,n_k,lj_k)      ! n=1,2,3,...
        rho_jl = OBD_Save(t_j,n_j,lj_j,t_l,n_l,lj_l)
        rho_il = OBD_Save(t_i,n_i,lj_i,t_l,n_l,lj_l)
        rho_jk = OBD_Save(t_j,n_j,lj_j,t_k,n_k,lj_k)
!     ...................................................................
        iphase = (-1)**(JJ+(j_k+j_l)/2+1)
        dme  = zzero
        d_ik = 0
        d_il = 0
        d_jk = 0
        d_jl = 0
        if(n_i.eq.n_k.and.t_i.eq.t_k.and.l_i.eq.l_k.and.j_i.eq.j_k) d_ik = 1
        if(n_i.eq.n_l.and.t_i.eq.t_l.and.l_i.eq.l_l.and.j_i.eq.j_l) d_il = 1
        if(n_j.eq.n_k.and.t_j.eq.t_k.and.l_j.eq.l_k.and.j_j.eq.j_k) d_jk = 1
        if(n_j.eq.n_l.and.t_j.eq.t_l.and.l_j.eq.l_l.and.j_j.eq.j_l) d_jl = 1

        dme   =              d_ik*d_jl*rho_ik*rho_jl           &
     &           +    iphase*d_il*d_jk*rho_il*rho_jk


!  -------------------- two-body density matrix elemetns of core states
!    For valence states
!  ........................................
         if( SPB%VType(lj_i,n_i-1,t_i).eq.1.and. SPB%VType(lj_j,n_j-1,t_j).eq.1    &
     & .and. SPB%VType(lj_k,n_k-1,t_k).eq.1.and. SPB%VType(lj_l,n_l-1,t_l).eq.1) then

             TBD   = Dens%ZRho2BJJ(a_ij,a_kl,bb) !+ ZRho2BJJ(a_kl,a_ij,bb))/2.d0  
        else
!  ........................................
!    For core nucleon  states and coupling between core and valence states 
!  ........................................
             TBD = dme*znorm                 ! (ij)<= (kl) 

        endif
!  ..........................................................
         Rho2BN = zzero
         if(abs(znorm) .gt. 1.d-100)  Rho2BN = TBD/znorm
!

          if(mod(tt,2) .eq.0 .and. a_ij.eq.a_kl) &
     &    Trace(tt/2) = Trace(tt/2) + (2.d0*JJ+1.d0)*Rho2BN*Metric

          if(abs(Rho2BN).gt.CHOP) write(ie,198) b_ij,b_kl,dreal(Rho2BN) 



           if(abs(Rho2BN).gt.1.d-8) write(300,'(3i6,f20.10)')        &
     &     a_ij,a_kl,bb,dreal(Rho2BN) 

         enddo
         enddo
            write(ie,20)
         enddo
!      .....................
         close(ie)
         close(300)
         write(*,'(a,i2)') '....The 2B-DME has been written to fort.',ie
         write(*,'(a,i2)') '....The 2B-DME is alo written to Rho2B.dat'

       write(*,*) 'Trace[Rho2B](n)=',Trace(0)
       write(*,*) 'Trace[Rho2B](p)=',Trace(1)
198    format(2i10,f20.10)
       end



        subroutine ZRho2B2Lambda2B(ie,ZRho1BJ0,ZLambda2B)
        USE VAPHFB_PAR
        implicit none
        integer   ie
        integer   J,k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        integer   bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        complex*16 zff,ZLambda2B,znorm,Rho2BN,Trace(0:1)
        Complex*16  ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        !DIMENSION ZLambda2B  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        DIMENSION ZLambda2B  (0:aMaxMax,0:aMaxMax,0:bMax)
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
        real*8  Metric
        integer  d_12,d_34,iphase
        integer  d_ik,d_jl,d_il,d_jk
        complex*16 rho_ik,rho_jl,rho_il,rho_jk
        complex*16 dme,TBD
!     ...............................

      ZLambda2B  (:,:,:) =zzero
      open(ie,file='Lambda2B.dat',status='unknown')
 
 20   format(/)
 50   format('*** Lamdab2B from GCM calculation ***')
 60   format('# basis (block, aMax)')
195   format('bMax=',i6)
196   format(2i10)
197   format('# block',i6,':',' aMax=',i6)
      write(ie,50)
      write(ie,60)
      write(ie,195) TPB%bMax-1
      if(TPB%bMax.gt.bMax) stop 'bMax is too small'
      do bm=0, TPB%bMax-1          ! The block starts from 0 to bMax_calc-1 
         write(ie,196) bm,TPB%cMax(bm) ! right 
      enddo

      write(ie,20)

       Trace(0:1) = zzero
       do bm=0, TPB%bMax-1
           write(ie,197) bm,TPB%cMax(bm)
           tt= TPB%btt(bm)
           PP= TPB%bPP(bm)
           JJ= TPB%bJJ(bm)
           bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        do b_ij=0, TPB%aMax(bb)-1
           a_ij= b_ij +1
        t_i = TPB%t1(bb,a_ij)
        n_i = TPB%n1(bb,a_ij)  ! 
        l_i = TPB%twol1(bb,a_ij)  ! l_i is doubled
        j_i = TPB%j1(bb,a_ij)  ! j_i is doubled
        LJ_i = (l_i+j_i-1)/2

        t_j = TPB%t2(bb,a_ij)
        n_j = TPB%n2(bb,a_ij)  ! 
        l_j = TPB%twol2(bb,a_ij)  ! l_j is doubled
        j_j = TPB%j2(bb,a_ij)  ! j_j is doubled
        LJ_j =(l_j+j_j-1)/2

       if(a_ij .ne. TPB%a(bb,FSPB%tnlj(lj_i,n_i,t_i),FSPB%tnlj(lj_j,n_j,t_j))) &
     &     stop 'Error in QNs for TPB%a'

        do b_kl=0, TPB%aMax(bb)-1
           a_kl= b_kl  +1
        t_k = TPB%t1(bb,a_kl)
        n_k = TPB%n1(bb,a_kl)  ! 
        l_k = TPB%twol1(bb,a_kl)  ! l_i is doubled
        j_k = TPB%j1(bb,a_kl)  ! j_i is doubled
        LJ_k = (l_k+j_k-1)/2

        t_l = TPB%t2(bb,a_kl)
        n_l = TPB%n2(bb,a_kl)  ! 
        l_l = TPB%twol2(bb,a_kl)  ! l_j is doubled
        j_l = TPB%j2(bb,a_kl)  ! j_j is doubled
        LJ_l = (l_l+j_l-1)/2

       if(t_i.eq.t_j.and.n_i.eq.n_j.and.LJ_i.eq.LJ_j) then
             Metric=(1+(-1)**(JJ))/2
        else
             Metric=2
        endif

!     ............................................................
        rho_ik = zzero
        rho_il = zzero
        rho_jk = zzero
        rho_jl = zzero

!        if(t_i.eq.t_k .and. lj_i.eq.lj_k) Rho_ik = ZRho1BJ0(t_i,lj_i,n_i-1,n_k-1)    ! n=0,1,2
!        if(t_i.eq.t_l .and. lj_i.eq.lj_l) Rho_il = ZRho1BJ0(t_i,lj_i,n_i-1,n_l-1)    ! n=0,1,2
!        if(t_j.eq.t_k .and. lj_j.eq.lj_k) Rho_jk = ZRho1BJ0(t_j,lj_j,n_j-1,n_k-1)    ! n=0,1,2
!        if(t_j.eq.t_l .and. lj_j.eq.lj_l) Rho_jl = ZRho1BJ0(t_j,lj_j,n_j-1,n_l-1)    ! n=0,1,2

!     ...................................................................
        iphase = (-1)**(JJ+(j_k+j_l)/2+1)
        dme  = zzero
        d_ik = 0
        d_il = 0
        d_jk = 0
        d_jl = 0
        if(n_i.eq.n_k.and.t_i.eq.t_k.and.l_i.eq.l_k.and.j_i.eq.j_k) d_ik = 1
        if(n_i.eq.n_l.and.t_i.eq.t_l.and.l_i.eq.l_l.and.j_i.eq.j_l) d_il = 1
        if(n_j.eq.n_k.and.t_j.eq.t_k.and.l_j.eq.l_k.and.j_j.eq.j_k) d_jk = 1
        if(n_j.eq.n_l.and.t_j.eq.t_l.and.l_j.eq.l_l.and.j_j.eq.j_l) d_jl = 1

!        dme   =              d_ik*d_jl*rho_ik*rho_jl           &
!     &           +    iphase*d_il*d_jk*rho_il*rho_jk

        dme   = d_ik*d_jl*ZRho1BJ0(t_i,lj_i,n_i-1,n_k-1)*ZRho1BJ0(t_j,lj_j,n_j-1,n_l-1)           &
     &        + iphase*d_il*d_jk*ZRho1BJ0(t_i,lj_i,n_i-1,n_l-1)*ZRho1BJ0(t_j,lj_j,n_j-1,n_k-1)

          Rho2BN   = Dens%ZRho2BJJ(a_ij,a_kl,bb) - dme   
          ZLambda2B(a_ij,a_kl,bb) = Rho2BN
!  .....................................................
!

          if(mod(tt,2) .eq.0 .and. a_ij.eq.a_kl) &
     &    Trace(tt/2) = Trace(tt/2) + (2.d0*JJ+1.d0)*Rho2BN*Metric

          if(abs(Rho2BN).gt.CHOP) write(ie,198)    b_ij,b_kl,dreal(Rho2BN) 

          if(bm.eq.0 .and. JJ.eq.0 .and. b_ij.eq.13 .and. b_kl.eq.43) then
           write(*,'(3i5)')  n_i,l_i,j_i
           write(*,'(3i5)')  n_j,l_j,j_j
           write(*,'(3i5)')  n_k,l_k,j_k
           write(*,'(3i5)')  n_l,l_l,j_l
          endif

         enddo
         enddo
            write(ie,20)
         enddo
!      .....................
         close(ie)
         write(*,'(a,i2)') 'The L2B has been written to Lambda2B.dat'

       write(*,*) 'Trace[L2B](n)=',Trace(0)
       write(*,*) 'Trace[L2B](p)=',Trace(1)
198    format(2i10,f20.10)
       end


!        subroutine Read_ME2B(ie,ZRho2BJJ)
        subroutine Read_ME2B(ie)
        USE VAPHFB_PAR
        implicit none
        Complex*16  me
        integer     ie,a12,a34,bb
!        Complex*16  ZRho2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)

        open(ie,file='Rho2B4Lambda2B.dat',status='old')
        read(ie,*)
        read(ie,*)
11      read(ie,*,end=21) a12,a34,bb,me

!        if(abs(me-2.0).lt.0.1) me = ztwo
!        if(abs(me-one).lt.0.1) me = zone
!        if(abs(me).lt.0.1)     me = zzero 
        Dens%ZRho2BJJ(a12,a34,bb) = me
        goto 11

21      continue
        call Check_Trace_ZRho2B()
        end


        subroutine Check_Trace_ZRho2B()
        USE VAPHFB_PAR
        implicit none
        integer   J,k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        integer   bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        complex*16 Trace(0:1)
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
        real*8  Metric
        INTEGER nn,np
!     ...............................
       Trace(0:1) = zzero
       do bm=0, TPB%bMax-1
           tt= TPB%btt(bm)
           PP= TPB%bPP(bm)
           JJ= TPB%bJJ(bm)
           bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        do b_ij=0, TPB%aMax(bb)-1
           a_ij= b_ij +1
        t_i = TPB%t1(bb,a_ij)
        n_i = TPB%n1(bb,a_ij)  ! 
        l_i = TPB%twol1(bb,a_ij)  ! l_i is doubled
        j_i = TPB%j1(bb,a_ij)  ! j_i is doubled
        LJ_i = (l_i+j_i-1)/2

        t_j = TPB%t2(bb,a_ij)
        n_j = TPB%n2(bb,a_ij)  ! 
        l_j = TPB%twol2(bb,a_ij)  ! l_j is doubled
        j_j = TPB%j2(bb,a_ij)  ! j_j is doubled
        LJ_j =(l_j+j_j-1)/2

           b_kl= b_ij
           a_kl= b_kl  +1
        t_k = TPB%t1(bb,a_kl)
        n_k = TPB%n1(bb,a_kl)  ! 
        l_k = TPB%twol1(bb,a_kl)  ! l_i is doubled
        j_k = TPB%j1(bb,a_kl)  ! j_i is doubled
        LJ_k = (l_k+j_k-1)/2

        t_l = TPB%t2(bb,a_kl)
        n_l = TPB%n2(bb,a_kl)  ! 
        l_l = TPB%twol2(bb,a_kl)  ! l_j is doubled
        j_l = TPB%j2(bb,a_kl)  ! j_j is doubled
        LJ_l = (l_l+j_l-1)/2

       if(t_i.eq.t_j.and.n_i.eq.n_j.and.LJ_i.eq.LJ_j) then
             Metric=(1+(-1)**(JJ))/2
        else
             Metric=2
        endif


          if(mod(tt,2) .eq.0 ) &
     &    Trace(tt/2) = Trace(tt/2) &
     &               + (2.d0*JJ+1.d0)*Metric*Dens%ZRho2BJJ(a_ij,a_kl,bb)

         enddo
         enddo

        nn = Nucl%nucleon(0) - Nucl%ncore(0)
        np = Nucl%nucleon(1) - Nucl%ncore(1)
        if(abs(Trace(0)-nn*(nn-1)).lt.1.d-3) then
           write(*,*) 'Trace[Rho2B](nn) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho2B](nn)=',Trace(0)
        endif
        if(abs(Trace(1)-np*(np-1)).lt.1.d-3) then
           write(*,*) 'Trace[Rho2B](pp) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho2B](pp)=',Trace(1)
        endif

       end
