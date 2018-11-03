       subroutine Iteration(NLEV,NCONS,AZ_P,AN_P,conv)
!      ...........................................................
!      step/velocity v: Z20 = - eta * (H_20 + Q_20 + ... )
!      ...........................................................
       use VAPHFB_Par
       implicit none
       integer III,ii,jj,imm,mm
       integer NLEV,NCONS
       real*8 AJZ2_MV,H_0,AZ_0,AN_0,AZ_P,AN_P,graddd
       real*8 Z20New(1:HO%NLEV*HO%NLEV),Z20(1:HO%NLEV*HO%NLEV)
       real*8 H20_0(1:HO%NLEV*HO%NLEV), O20(1:HO%NLEV*HO%NLEV)
       real*8 StNew,St,Z20sq,znorm
       real*8 MtNew(1:HO%NLEV*HO%NLEV),Mt(1:HO%NLEV*HO%NLEV)
       real*8, DIMENSION(:), allocatable :: ACONS_20 !(NCONS*NLEV*NLEV)
       real*8 bet2t,gam2t,P00t
       real*8 cj,cjxp,cjxn,omegax,conv,zg_0,EPNP
!      ...........................................................
       print *, ' -> Iteration ...'
       print *, '    ......... '
       write(*,'(A4,100A14)') 'It','Grad','E_HFB', 'Prot',&
       &'Neut','bet2','gam2','P00','   <J>x','   <J^2_z>'

        if(NCONS.gt.NCONSMAX) &
       & stop 'Error: NCONS is larger than NCONSMAX!'

!      ......... allocate memory
        if(.NOT. ALLOCATED(HFB%RO_0))      &
       &  ALLOCATE(HFB%RO_0(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%Akapa10_0)) &
       & ALLOCATE(HFB%Akapa10_0(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(HFB%Akapa01_0)) &
       & ALLOCATE(HFB%Akapa01_0(1:HO%NLEV,1:HO%NLEV))
        if(.NOT. ALLOCATED(ACONS_20))    &
       & ALLOCATE(ACONS_20(1:HO%NLEV*HO%NLEV*NCONS))

!        write(*,'(a,f6.3,a)')' Iteration: allocated ACONS_20 &
!     &  with memory:',sizeof(ACONS_20)/(1024*1024*1024.0),'G byte ..'
       

        ! initialization
        Z20sq = zero 
        znorm = 1.d-3
        St    = 1.d0 
        StNew = 1.d0 
        do jj=1,NLEV**2
           Z20(jj)=zero
           Mt(jj)=zero
        end do
!       ........................ start the iteration
        DO III=1,Input%itermax
!       .........................Initialitation
           H_0 =zero
           AN_0=zero
           AZ_0=zero

           do jj=1,NLEV**2
              H20_0(jj)=zero
           end do
!          write(*,*) ' ..... Rho and kappa ...'
           call UV2Density(HFB%U0,HFB%V0,HFB%RO_0,HFB%Akapa10_0,HFB%Akapa01_0,HO%NLEV)
           call N_0(HFB%V0,AN_0,AZ_0,HO%NLEV)
           call betgamm_0(HFB%V0,bet2t,gam2t,HO%NLEV)
           call Qpair_0(cME1B%P_00_10_me,HFB%Akapa10_0,HFB%Akapa01_0,P00t,HO%NLEV)

           do mm=1,NCONS-ipair_const
              imm=(mm-1)*HO%NLEV**2+1
              call F_20(HFB%U0,HFB%V0,ACONS_ME(imm),ACONS_20(imm),NLEV) 
           end do

           do mm=NCONS-ipair_const+1,NCONS
              imm=(mm-1)*HO%NLEV**2+1
              call G_20(HFB%U0,HFB%V0,ACONS_ME(imm),ACONS_20(imm),NLEV)
           end do
!       ........................................
           call NZPROJ_VAP(AZ_P,AN_P,EPNP,H20_0,PNP%NFOM,NLEV)

        if (III.eq.1) then
           do mm=1,NCONS
              alag_0(mm)=0.d0
           end do
        end if

!      .....

!        write(*,*) '.... LAGRANGE_MULT ........'
        call LAGRANGE_MULT(NCONS,H20_0,ACONS_20,alag_0,alag,NLEV)

        HFB%EFermi(0) =  alag_0(1)
        HFB%EFermi(1) =  alag_0(2)
        do mm=1,NCONS
           write(111,*) iii, alag(mm)
        enddo
!.......Gradient method.......
        Z20sq = zero        
        do jj=1,NLEV**2
          graddd=0.d0
          do mm=1,NCONS
             imm=(mm-1)*NLEV**2
             graddd=graddd+alag(mm)*ACONS_20(imm+jj)
          end do
          O20(jj) = H20_0(jj)-graddd 
          Z20sq   = Z20sq + (H20_0(jj)-graddd)**2        


        ! optimizers
        ! https://distill.pub/2017/momentum/
!        HFB%Thoul_0(jj)=  -Input%eta*(H20_0(jj)-graddd) !aqui
        if(Opt%method == 'GD') then
           Z20New(jj)=  Opt%lr*(H20_0(jj)-graddd) !aqui
        else if(Opt%method == 'GDM') then
           Z20New(jj)=  Opt%gam1*Z20(jj) + Opt%lr*O20(jj)
        else if(Opt%method == 'ADAM') then
        else
           print *, ' Optimizer is NOT implemented '
        endif
        HFB%Thoul_0(jj) = -Z20New(jj)
        end do

        if(Opt%method == 'ADAM') then
        do jj=1,NLEV**2
          StNew     = Opt%gam2*St + Z20sq*1.d-4
          !StNew     = Opt%bet2*St/(1.d0-Opt%bet2) + Z20sq*znorm*0.1
          MtNew(jj) =  Opt%gam1*Mt(jj) + O20(jj)
          Z20New(jj)=  MtNew(jj)*Opt%lr/(sqrt(Stnew)+Opt%eps)
          HFB%Thoul_0(jj) = -Z20New(jj)
        end do
        !print *, ' sq(Z20sq*znorm) =',sqrt(Z20sq*znorm)
        !print *, ' adapted step size =',Opt%lr/(sqrt(Stnew)+Opt%eps)
        endif

        !  C := alpha*op( A )*op( B ) + beta*C,
        call DGEMM ('t','n',NLEV,NLEV,NLEV,1.d0,HFB%Thoul_0,NLEV,HFB%Thoul_0,NLEV,0.d0,HFB%GRAD_0,NLEV)

        zg_0=zero
        do ii=1,NLEV
           zg_0=zg_0+HFB%GRAD_0(ii,ii)
        end do

        conv=zg_0*StNew/(Opt%lr**2)

!      .... check <J_x>
        call Q_0(cME1B%AJX_ME,HFB%RO_0,cjxp,cjxn,NLEV)
!      .... check <J2_z>
        call F2_0(HFB%RO_0,HFB%AKAPA10_0,HFB%AKAPA01_0,  &
     &            cME1B%AJZ_ME,cME1B%AJZ2_ME,AJZ2_MV,NLEV)
        cj = (sqrt(1.0+4.0*(cjxp+cjxn)**2)-1)*0.5

        write(*,1001),III,conv,EPNP+H%E0,AZ_P,AN_P,bet2t,gam2t,P00t,cj,AJZ2_MV
        if (conv.le.Input%tolgrad) then
           if(Input%bkick.lt.CHOP) then
             print*,'CALCULATION CONVERGED'
             print*,' '
             exit
!      ..........................................
            elseif(Input%bkick.gt.CHOP .and. HFB%iconv.gt.0) then  ! with a nonzero kick stepsize 
              print*,'CALCULATION CONVERGED'
              print*,' '
              exit
!     .................... after the first convergence     
            elseif(Input%bkick.gt.CHOP .and. HFB%iconv.eq.0) then  ! with a nonzero kick stepsize 
              HFB%Et = EPNP+H%E0    ! record the energy
              print*,' kick the state with a Delta_bet2=',Input%bkick
!     .................... record the U,V      
              do ii=1,NLEV
              do jj=1,NLEV
                 HFB%Ut(jj,ii)=HFB%U0(jj,ii)
                 HFB%Vt(jj,ii)=HFB%V0(jj,ii)
              end do
              end do
! add a kick onto Q20 after the first convergence
              if(Const%iQB.eq.1) acons_mv(3) = qtc20+Input%bkick/Const%Q2BA
              if(Const%iQB.ne.1) acons_mv(3) = qtc20+Input%bkick
              HFB%iconv   = HFB%iconv + 1
           endif
        endif
1001    format(I4,100F14.6)


        call adjust_cons(NCONS,HFB%U0,HFB%V0,HFB%Thoul_0,HFB%U1,HFB%V1,HFB%ALcholes,NLEV)


!      update U,V       
        do ii=1,NLEV
         do jj=1,NLEV
            HFB%U0(jj,ii)=HFB%U1(jj,ii)
            HFB%V0(jj,ii)=HFB%V1(jj,ii)
         end do
        end do

!      update quantities for optimiziers
        St = StNew       
        do jj=1,NLEV**2
!            HFB%Thoul_old(jj)=HFB%Thoul_0(jj)
            ! store the velocity/gradient
            Z20(jj) = Z20New(jj)
            Mt(jj)  = MtNew(jj)
!            HFB%Thoul_0(jj)  =zero
        end do
        if(III.eq.1) znorm = 1.d0/Z20sq

        !stop 'I was stop in lib/HFBA_Iteration !'
!      ........................ print out intermediate wave functions
        if(mod(III,20).eq.0) then
        !writing the wave function of check point
         print *, 'write check point file ..'
        write(12,*) HFB%ide
         do ii=1,NLEV
          do jj=1,NLEV
             write(12,*) HFB%U0(ii,jj),HFB%V0(ii,jj)
          end do
         end do
         close(12)
        endif

        if(III.eq.Input%itermax) then
         print*,'MAXIMUM NUMBER OF ITERATIONS REACHED'
         print*,' '
        end if

!       .................................................
        END DO
        !END OF THE GRADIENT ITERATIONS
!       .................................................
        if(Input%bkick.gt.CHOP .and. HFB%Et .lt. EPNP+H%E0) then ! keep the U,V of the first convergence  
        print *,' recover the wave function from the first convergence.'
        do ii=1,NLEV
         do jj=1,NLEV
            HFB%U0(jj,ii)=HFB%Ut(jj,ii)
            HFB%V0(jj,ii)=HFB%Vt(jj,ii)
         end do
        end do
        endif

       return
       end
