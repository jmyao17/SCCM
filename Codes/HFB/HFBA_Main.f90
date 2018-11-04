!      .....................................................
       program HFB_MPI
!      .....................................................
!      HFB solver:
!      1) Particle Number Projection Before Variation
!      2) Constraints: 
!         quadrupole deformations, 
!         iso-scalar pairing,
!         iso-vector pairing,
!         etc.
!      3) A kick along Q20 is implemented after the 1st convergence 
!         to jump out from the possible local minimum
!      4) pairing collaps sometimes only for either neutrons or protons
!         (not for both). In this case, it is better to constrain a little bit
!         (<P0>=0.2) np isoscalar pairing amplitude such that the Pfaffian method
!         for the norm overlap in the PNAMP part will still work.  
!      ....................................................
       use mpi

       USE VAPHFB_Par
       implicit none
       integer ii,jj,iq,iqout
       real*8  AZ_P,AN_P,epsi

       integer :: ierr,nprocs,rank  ! variables for MPI
       logical :: lexist
       !
       ! initialize mpi
       call MPI_INIT(ierr)
       ! setup communicator size
       call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
       ! setup ranks/IDs for each process
       call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)


       if(rank.eq.0) then
         open(82,file='E_beta_gam.dat',status='unknown')
         write(82,'(a2,9A8,6A10)') '','epsi','cf','bet2n','bet2p', 'bet2t',&
        'gam2n','gam2p','gam2t',"P_00",'EP_NN','EP_PP','EP_NP','Etot','Radius'
       endif

       if(rank.eq.0) call Reader(.true.)
       if(rank.ne.0) call Reader(.false.)
       if(rank.eq.0) call Prep(.true.)
       if(rank.ne.0) call Prep(.false.)
       call Model_space()
       call ZCONST
       call gfv
       call Hamiltonian(rank,.false.) 
       call Init_hfbwf()
       CALL ZNme()
       call Isospin_t3()


      if(rank .eq.0) print *, ' -> Constraint  '
      if(Const%NCONS.lt.2) stop ' Error: Const%NCONS value is wrong !'
      if(.NOT. ALLOCATED(ACONS_ME))  ALLOCATE(ACONS_ME(1:Const%NCONS*HO%NLEV**2))
      if(rank .eq.0)  write(*,'(a,f8.3,a)')'     ME1B4Constraint: allocated ACONS_ME &
     &with memory:', sizeof(ACONS_ME)/(1024*1024*1024.0),'G byte ..'

        do ii=1,HO%NLEV**2
           ACONS_ME(ii+0*HO%NLEV**2)=AZme(ii)
           ACONS_ME(ii+1*HO%NLEV**2)=ANme(ii)
        end do

        if(rank .eq.0) write(*,*) '    Prepare 1B constraint matrix elements O_ab ...'
       call ME1B4Constraint()
        if(rank .eq.0) write(*,*) '    Read constraint coordinates <O_i> ...'
       call ReadGeneCoord()

      if(.NOT. ALLOCATED(HFB%gamma_0))  ALLOCATE(HFB%gamma_0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%delta10))  ALLOCATE(HFB%delta10(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%ham_0))    ALLOCATE(HFB%ham_0(1:HO%NLEV,1:HO%NLEV))



!     .........................................
!     Loop over coordinates: beta,gamma,p00
!     .........................................
        if(rank .eq.0) write(*,*) '    Loops over the constraint coordinates <O_i> ...'
      do iq=1,Const%NQ

         if(rank .eq. mod(iq,nprocs)) then
         !if(rank .eq. iq) then
        
         
         HFB%iconv=0
         acons_mv_aux(12)  = Const%P00_mesh(iq)
         shared%cf = Const%hw(iq)
         write(*,*) '    ..........................................'
         write(*,'(a,f6.3,a,f8.3,a,f6.3)') &
    &       '   beta=',Const%beta2t_mesh(iq),&
    &       '   gamma=',Const%gamma2t_mesh(iq), &
    &       '   P00=',Const%P00_mesh(iq), &
    &       '   crank_freq=',Const%hw(iq)
       write(*,*) '    ..........................................'
       call Filename4wfs(Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),&
                         Const%P00_mesh(iq),Const%hw(iq))
!      ........................ read wf
       if(Input%inwf.eq.1) then
         inquire (file=File%wf,exist=lexist)
         if(.not.lexist) goto 8
         open(10,file=File%wf,err=8)
         write(*,*)'    read w.f. from:',File%wf
         read (10,*,err=8) HFB%ide
         do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             read(10,*,end=8) HFB%U0(ii,jj),HFB%V0(ii,jj)
          end do
         end do
  8      continue
         close(10)
       endif

 9    print *,'   Set Constraint terms ...' 
       call SetConstraintValue(iq,HO%NLEV)

      if(.NOT. ALLOCATED(HFB%Thoul_0_work))  ALLOCATE(HFB%Thoul_0_work(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%Thoul_0))       ALLOCATE(HFB%Thoul_0(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(HFB%Thoul_old))     ALLOCATE(HFB%Thoul_old(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(HFB%GRAD_0))  ALLOCATE(HFB%GRAD_0(1:HO%NLEV,1:HO%NLEV))
        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            HFB%Thoul_0_work(ii,jj)=zero
         end do
        end do

!        write(*,*) 'Main: adjust_cons'
        if(.NOT. ALLOCATED(HFB%U1))  ALLOCATE(HFB%U1(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%V1))  ALLOCATE(HFB%V1(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%ALcholes))  ALLOCATE(HFB%ALcholes(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%Ut))  ALLOCATE(HFB%Ut(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%Vt))  ALLOCATE(HFB%Vt(1:HO%NLEV,1:HO%NLEV))

        print *,'   adjust constraint terms ...' 
        call adjust_cons(Const%Ncons,HFB%U0,HFB%V0,HFB%Thoul_0_work,HFB%U1,HFB%V1,&
      &                  HFB%ALcholes,HO%NLEV)


        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
          HFB%U0(jj,ii)=HFB%U1(jj,ii)
          HFB%V0(jj,ii)=HFB%V1(jj,ii)
          HFB%Thoul_0_work(ii,jj)=zero
!          write(231,*) HFB%U0(jj,ii),HFB%V0(jj,ii)
         end do
        end do
!       .................................................
!       gradient method is adopted

        call Iteration(HO%NLEV,Const%Ncons,AZ_P,AN_P,epsi) 
!       .................................................

  
        !writing the wave function
!         write(11,*) HFB%ide
!         do ii=1,HO%NLEV
!          do jj=1,HO%NLEV
!           write(11,*) HFB%U0(ii,jj),HFB%V0(ii,jj)
!          end do
!         end do
!         print *, ' .. checking U, V ...'
!         call checkwf(HFB%U0,HFB%V0,HO%NLEV)
!       .................................................
        !OUTPUT
!       .................................................
        iqout = iq*5
       open(iqout,file=File%out,status='unknown')
       call Calc_Obs_From_UV(rank,iqout,epsi,HFB%U0,HFB%V0,HFB%ro_0,HFB%akapa10_0,HFB%akapa01_0, &
     &             H%ME1BM,AZ_P,AN_P,cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t,    &
     &             cME1B%AJX_ME,cME1B%AJX2_ME,cME1B%AJY_ME,cME1B%AJY2_ME,cME1B%AJZ_ME,cME1B%AJZ2_ME, &
     &             cME1B%P_1m1_00_me,cME1B%P_10_00_me,cME1B%P_1p1_00_me, &
     &             cME1B%P_00_1m1_me,cME1B%P_00_10_me,cME1B%P_00_1p1_me,HO%NLEV,HO%NLEV**2)

       print *, ' -> Canonic basis  '
       call canon_basis(iqout,HFB%ro_0,HFB%gamma_0,HFB%delta10,HFB%ham_0,HO%NLEV)

       !call HFB_Diag(HFB%ham_0,HFB%delta10,HO%NLEV)

        open(11,file=File%wf,status='unknown')
        print *,' writting out wave function to ',File%wf
        if( HFB%ide(0)+ HFB%ide(1).eq.1) write(120,*) 'WARNING: ',File%wf, HFB%ide(0),HFB%ide(1)    
        write(*,*) HFB%ide
        write(11,*) HFB%ide
         do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             write(11,*) HFB%U0(ii,jj),HFB%V0(ii,jj)
          end do
         end do
       close(iqout)
!    .......................... loop over collective coordinate
       endif

       enddo ! iq

       print *, 'Rank:',rank,' is finished ...'       
       call MPI_FINALIZE(ierr)

       end
