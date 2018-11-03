        subroutine Hamiltonian(rank,lpr) 
!     ..........................
        USE VAPHFB_PAR
        use omp_lib
        implicit none
        logical lpr
        integer rank
        integer mitf,nn
        real*8, dimension(:,:,:,:,:,:,:), allocatable :: V_AB_CD_JTMT
        real*8, dimension(:,:,:,:,:,:), allocatable :: V_AB_CD_J,V_AB_CD_JT
        real*8 start,finish
        integer ii,kk,mj,kk_n
        integer ia,ib,ic,id
        integer mjinc,jttot(1:2)
        integer lang(HO%NLEV),jang(HO%NLEV),mjang(HO%NLEV),mtisos(HO%NLEV) !uncoupled basis
        integer jindex(HO%NLEV),nlindex(HO%NLEV)
        integer jmax,nlj,iabcd
        real*8  hja,hjb,hjc,hjd
        integer it12,MMT,iJ,iT
        real*8  ajtot,attot,ame_coup
        real*8  phasJT,phasab,phascd,delta_ab,delta_cd
        integer jja,jma,jta,jna,jlja,jacoup
        integer jjb,jmb,jtb,jnb,jljb,jbcoup
        integer jjc,jmc,jtc,jnc,jljc,jccoup
        integer jjd,jmd,jtd,jnd,jljd,jdcoup
        integer JJ,MMJ,MT
        integer icheck
        real*8  AN_INV,AN_AB,AN_CD,cb1,cb2,cb3,cb4,suma,sumando
        real*8  check,cg1
        character(500) find_file
        character*14 ctemp
!      .................... precalculate the CG-coefficient
        call PreCalc_CG()

       if(Input%IntIMSRG .ne. 1) then
          stop 'IntIMSRG should be 1'
       endif

!       if(Input%IntType.ne.1) then
!          File%f1b="../Int/IMSRG_"//trim(Input%cFlow)//"_SPE_"//trim(Input%cIntID)//&
!          & "_"//trim(Input%cValID)//".dat"
!       endif
!       if(Input%IntType.eq.1) then
          File%f1b="../Int/IMSRG_"//trim(Input%cFlow)//"_SPE_"//trim(Input%cIntID)//&
          & "_"//trim(Input%cValID)//'_'//Input%chwHO//".dat"
!       endif
!      ..... two-body
       if(Input%IntType.eq.1) then
           File%f2bJ='../Int/IMSRG_'//trim(Input%cFlow)// &
        &  '_'//trim(Input%cIntID)//trim(Input%ctpp)//trim(Input%cValID) &
        &  //'_'//Input%chwHO//'_J.dat'
          else
           File%f2bJ='../Int/IMSRG_'//trim(Input%cFlow)// &
        &  '_'//trim(Input%cIntID)//'_'//trim(Input%cValID) &
        &  //'_'//Input%chwHO//'_J.dat'
       endif         

       if(Input%IntType.eq.1 ) then
           File%f2bJTMT='../Int/IMSRG_'//trim(Input%cFlow)//'_'//trim(Input%cIntID)//&
        &  Input%ctpp//trim(Input%cValID)//'_'//Input%chwHO//'_JTMT.dat'
        endif
       if(Input%IntType.eq.0) then
          File%f2bJTMT='../Int/IMSRG_'//trim(Input%cFlow)//'_'&
     &    //trim(Input%cIntID)//'_'//trim(Input%cValID)//'_JTMT.dat'
       endif

       if(Input%IntType.eq.1) then
         File%f2bm='IMSRG_'//trim(Input%cFlow)//'_'//trim(Input%cIntID)//&
         &Input%ctpp//trim(Input%cValID)//'_'//Input%chwHO//&
         &'_m1m2m3m4.dat'
       else
         File%f2bm='IMSRG_'//trim(Input%cFlow)//'_'//trim(Input%cIntID)//&
         &'_m1m2m3m4.dat'
       endif
       print *,'rank:',rank,'->',Input%chwHO,File%f2bm !,INT_DIR
!      ....... remove spaces in the filenames
      
       call StripSpaces(File%f1b)
       call StripSpaces(File%f2bm)
       call StripSpaces(File%f2bJTMT)

       if(Input%IntJT.eq.0) INT_DIR = find_file("GCM_ME_FILES",trim(File%f2bm))
      
!       if(Input%IntJT .ne.0) then
          open(41,file='chi2b3b.int',status='old')
          read(41,'(a14,a)') ctemp,File%IMSRG_Hme1b
          read(41,'(a14,a)') ctemp,File%IMSRG_Hme2b
          close(41)
          if(lpr) print *,File%IMSRG_Hme1b,File%IMSRG_Hme2b
          INT_DIR = find_file("GCM_ME_FILES",File%IMSRG_Hme1b)
!       endif

       if(lpr) then
         print *, ' ======================== '
         print *, ' Details of Calculations  '
         print *, ' ======================== '
         write(*,*) ' -----------------------'
         write(*,*) ' idx   n+1    l    2j '
         do nlj=1,HO%nljmax
!           nlshort(nlj)=SPB%n(nlj)+1   ! 1d5/2; n =1,2,..
            write(*,'(4i5)') nlj,SPB%n(nlj)+1,SPB%l(nlj),SPB%twoj(nlj)
          enddo    
        endif

      !possible values of T=t1+t2 
      jttot(1)=0
      jttot(2)=1


      !setting up the values |j m_j> 
      kk=0
      do ii=1,HO%nljmax       ! max number of nlj value
        jmax=SPB%twoj(ii)+1  ! 2j+1
        mjinc=0

        do mj=1,jmax      !jmax= 2j+1
          kk=kk+1
!       ......... proton
          nlindex(kk)=SPB%n(ii)+1 !nlshort(ii)
          lang(kk)   =SPB%l(ii) ! lshort(ii)
          jang(kk)   =SPB%twoj(ii) !jshort(ii)
          mjang(kk)  =SPB%twoj(ii) - mjinc
          mtisos(kk) =-1              ! isospin -1 for proton
          jindex(kk) = ii 
!       ......... neutron
          kk_n = kk + HO%NOrbit(0)
          nlindex(kk_n)=SPB%n(ii)+1 ! nlshort(ii)
          lang(kk_n)   =SPB%l(ii) ! lshort(ii)
          jang(kk_n)   =SPB%twoj(ii) !jshort(ii)
          mjang(kk_n)  =SPB%twoj(ii)-mjinc
          mtisos(kk_n) =1            ! isospin +1 for neutron
          jindex(kk_n) = ii 

          mjinc=mjinc+2
        end do
      end do
      if(kk_n .ne. HO%NLEV) then
         print *, 'HO%NLEV should be',kk_n
         stop 'Main: NLEV is wrong ...'
      endif
!     ............... quantum numbers for each orbital
      if(.NOT. ALLOCATED(tnljm%t)) ALLOCATE(tnljm%t(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%n)) ALLOCATE(tnljm%n(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%l)) ALLOCATE(tnljm%l(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%lj)) ALLOCATE(tnljm%lj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twoj)) ALLOCATE(tnljm%twoj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twom)) ALLOCATE(tnljm%twom(1:HO%NLEV))
      if(.NOT. ALLOCATED(HO%tts)) ALLOCATE(HO%tts(1:HO%NLEV))
       if(.NOT. ALLOCATED(tnljm%level)) &
     & ALLOCATE(tnljm%level(-HO%twojmax:HO%twojmax,0:HO%ljmax,0:HO%nmax,0:1))

!      if(lpr) then
!        write(*,*) '.................................'
!        write(*,*) '| Single-Particle Basis         |' 
!        write(*,*) '.................................'
!        write(*,*)  'k   t   n+1   l   twoj   twom'
!      endif
      do kk=1,kk_n       ! max number of s.p. state
          tnljm%t(kk) = mtisos(kk) 
          tnljm%n(kk) = nlindex(kk)
          tnljm%l(kk) = lang(kk)
          tnljm%twoj(kk) = jang(kk)
          tnljm%lj(kk)   = (2*lang(kk)+jang(kk)-1)/2
          tnljm%twom(kk) = mjang(kk)
         tnljm%level(mjang(kk),tnljm%lj(kk),tnljm%n(kk),mitf(mtisos(kk))) = kk
!       if(lpr) &   
!      &write(*,20)  kk,tnljm%t(kk),tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),'/2',tnljm%twom(kk),'/2'

! ............
          nn = 2*tnljm%n(kk)+tnljm%l(kk)
          if (nn.lt.10) then
              write(HO%tts(kk),19) tnljm%n(kk),HO%tl(tnljm%l(kk)), &
     &                 tnljm%twoj(kk),' /2 ',tnljm%twom(kk),' /2 ',    &
     &                 HO%tp(mip(iv(tnljm%l(kk)))) 
               else
                  write(HO%tts(kk),21) tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),tnljm%twom(kk)
               endif

      end do
19   format(1h[,i2,a1,i3,a4,i3,a3,1h],a2)
20   format(5i4,a,i4,a)
21   format(4i2)
!   .................................
      if(lpr) then  
        print *, ' ...............'
        print *, ' -> Hamiltonian '
      endif
!      one-body energy
       if(.NOT. ALLOCATED(H%ME1BM)) ALLOCATE(H%ME1BM(1:HO%NLEV,1:HO%NLEV))

!     ............... read ME1B from files 
      !if(Input%IntJT.eq.3) call Generate_ME1B(H%E0,H%ME1BM)   ! read ME1B from the files by Nathan
      !if(Input%IntJT.ne.3) 

       call kinetic(H%ME1BM,H%ddd_kin)


!      determine the dimension
 
      if(lpr) write(*,*) '........... determine the iabcd_max  ...'
      iabcd = 0
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(iv,HO,tnljm) REDUCTION(+: iabcd)
!$OMP DO SCHEDULE(DYNAMIC) 
      do ia=1,HO%NLEV
      do ib=ia+1,HO%NLEV
!     .................. do ib=ia+1,HO%NLEV
         do ic=1,HO%NLEV
!         ...............do id=1,HO%NLEV
         do id=ic+1,HO%NLEV
         if(tnljm%t(ia)+tnljm%t(ib) .ne. tnljm%t(ic)+tnljm%t(id))             cycle   ! isospin conservation
         if(tnljm%twom(ia)+tnljm%twom(ib) .ne. tnljm%twom(ic)+tnljm%twom(id)) cycle   ! m conservation
         if(iv(tnljm%l(ia)+tnljm%l(ib)) .ne. iv(tnljm%l(ic)+tnljm%l(id)))     cycle   ! parity
         if(ia+ib .gt. ic+id)                                                 cycle
         
         iabcd     = iabcd + 1
! ...................
      enddo
      enddo
      enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      H%iabcd_max = iabcd
     
!      print *, ' iabcd_max=',iabcd

      if(.NOT. ALLOCATED(H%ka))  ALLOCATE(H%ka(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kb))  ALLOCATE(H%kb(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kc))  ALLOCATE(H%kc(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kd))  ALLOCATE(H%kd(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%ME2BM))  ALLOCATE(H%ME2BM(1:H%iabcd_max))
      if(lpr) then
       write(*,*) '   No. of stored two-body matrix elements:',H%iabcd_max
       write(*,'(a,f6.3,a)')'    Main: allocated ka,kb,kc,kd with &
     & memory:', 4*sizeof(H%ka)/(1024*1024*1024.0),'G byte ..'
       write(*,'(a,f6.3,a)')'    Main: allocated H%ME2B with memory:', &
     & sizeof(H%ME2BM)/(1024*1024*1024.0),'G byte ..'
      endif

      if(lpr) print *, ' -> Initializing ka,kb,kc,kd ...'
      iabcd = 0
      icheck = HO%NLEV/100
      if(icheck.eq.0) icheck=HO%NLEV/10

!! Using OpenMP take even longer time 
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(iabcd,icheck,iv,HO,tnljm,H)
!!$OMP DO SCHEDULE(DYNAMIC)
      do ia=1,HO%NLEV
         if(lpr .and. mod(ia,icheck).eq.0) call progress_bar(int(ia/icheck))  
      do ib=ia+1,HO%NLEV
         do ic=1,HO%NLEV
         do id=ic+1,HO%NLEV
         if(tnljm%t(ia)+tnljm%t(ib) .ne. tnljm%t(ic)+tnljm%t(id))             cycle   ! isospin conservation
         if(tnljm%twom(ia)+tnljm%twom(ib) .ne. tnljm%twom(ic)+tnljm%twom(id)) cycle   ! m conservation
         if(iv(tnljm%l(ia)+tnljm%l(ib)) .ne. iv(tnljm%l(ic)+tnljm%l(id)))     cycle   ! parity
         if(ia+ib .gt. ic+id)                                                 cycle
!!$OMP CRITICAL
         iabcd     = iabcd + 1
         H%ka(iabcd) = ia
         H%kb(iabcd) = ib
         H%kc(iabcd) = ic
         H%kd(iabcd) = id
!         H%ME2BM(iabcd) = 0.d0
!!$OMP END CRITICAL
      enddo
      enddo
      enddo
      enddo
!!$OMP END DO
!!$OMP END PARALLEL

      H%ME2BM(:) = 0.d0
      if(lpr) then
        write(*,*)
        write(*,*) ' -> Read hamiltonian matrix elements'
      endif








!     ..................................................
!     decoupling ME2B from JTMT to m-scheme
!     ..................................................
      if(Input%IntJT.eq.3) then
!     ..................................................
!     ......... Initializing the two-body matrix elements in coupled basis
!       if(.NOT. ALLOCATED(V_AB_CD_J))  &
!     & ALLOCATE(V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5))
!       if(.NOT. ALLOCATED(V_AB_CD_JT))  &
!     & ALLOCATE(V_AB_CD_JT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1))
!      do ia=1,HO%nljmax
!       do ib=1,HO%nljmax
!        do ic=1,HO%nljmax
!         do id=1,HO%nljmax
!          do iJ=0,HO%twojmax
!  ...................
!           do iT=0,HO%tmax
!              V_AB_CD_JT(ia,ib,ic,id,iJ,iT) =0.d0
!           end do
!           do it12=0,5
!            V_AB_CD_J(ia,ib,ic,id,iJ,it12)=0.d0
!           end do
!  ...................
!          end do
!         end do
!        end do
!       end do
!      end do


       open(99,file=trim(INT_DIR)//trim(File%f2bm),&
               form='unformatted',status='unknown')

!  ...................................
!      READ H without Isospin Symmetry
!  ....................................

       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))

      !write(*,*) '   the 2B terms from ',trim(File%f2bJTMT)
       open(1,file=trim(File%f2bJTMT),status='old')
       !write(*,*) '   JT=1: the 2B terms from  fort.48'
       !open(1,file='fort.48',status='old')
6      read(1,*,end=7) ia,ib,ic,id,iJ,iT,MMT,ame_coup

!       write(200,'(6i5,f12.8)') ia,ib,ic,id,iJ-1,iT-1,ame_coup
       hja=SPB%twoj(ia)/2.d0
       hjb=SPB%twoj(ib)/2.d0
       hjc=SPB%twoj(ic)/2.d0
       hjd=SPB%twoj(id)/2.d0
       ajtot=1.d0*iJ ! jtot(iJ)
       attot=1.d0*iT ! jttot(iT)

       phasab=(-1.d0)**(hja+hjb+ajtot+attot)
       phascd=(-1.d0)**(hjc+hjd+ajtot+attot)

       V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,ic,id,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,id,ic,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ia,ib,id,ic,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
!
       V_AB_CD_JTMT(ic,id,ia,ib,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ia,ib,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ib,ia,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ic,id,ib,ia,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme

! !     end do
       goto 6
  7   continue
      close(1)

      !computing the two-body matrix elements in uncoupled basis

      if(lpr) &
      print *, ' uncoupling ME2B from J-scheme to m-scheme and store to',trim(INT_DIR)//trim(File%f2bm)
!     ..............................................................................
!      OMP
!     ..............................................................................
      icheck = H%iabcd_max/100

!$OMP PARALLEL DEFAULT(PRIVATE)
!SHARED(icheck,CG_Save,H,HO,tnljm,jindex,V_AB_CD_JTMT) 
!$OMP DO schedule(dynamic)
      do iabcd = 1, H%iabcd_max
        ia = H%ka(iabcd)
        jja=tnljm%twoj(ia)  ! jang(ia)
        jma=tnljm%twom(ia)
        jta=tnljm%t(ia)      !mtisos(ia) 
        jna=tnljm%n(ia)
        jlja=tnljm%lj(ia)
        jacoup=jindex(ia)

!        print*,ia
        ib = H%kb(iabcd)
        jjb=tnljm%twoj(ib) ! jang(ib)
        jmb=tnljm%twom(ib)  ! mjang(ib)
        jtb=tnljm%t(ib)  !mtisos(ib)
        jnb=tnljm%n(ib)
        jljb=tnljm%lj(ib)
        jbcoup=jindex(ib)
        delta_ab=0.d0
        if(jna.eq.jnb .and. jlja.eq.jljb) delta_ab=1.d0
!        do ic=1,NLEV
        ic = H%kc(iabcd)
         jjc=tnljm%twoj(ic)   ! jang(ic)
         jmc=tnljm%twom(ic)   ! mjang(ic)
         jtc=tnljm%t(ic)      ! mtisos(ic)
         jnc=tnljm%n(ic)
         jljc=tnljm%lj(ic)
         jccoup=jindex(ic)
          id = H%kd(iabcd)
          jjd=tnljm%twoj(id)   !jang(id)
          jmd=tnljm%twom(id)   !mjang(id)
          jtd=tnljm%t(id)      !mtisos(id)
          jnd=tnljm%n(id)
          jljd=tnljm%lj(id)
          jdcoup=jindex(id)
          delta_cd=0.d0
        if(jnc.eq.jnd .and. jljc.eq.jljd) delta_cd=1.d0

          suma=0.d0
          DO JJ=0,HO%twojmax
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ

!              call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)
!              call CJJ(jjc,jjd,2*JJ,jmc,jmd,2*MMJ,cb2)
              cb1  = CG_Save((jja+1)/2,(jma+1)/2,(jjb+1)/2,(jmb+1)/2,JJ,MMJ) ! J12,M12 are NOT doubled 
              cb2  = CG_Save((jjc+1)/2,(jmc+1)/2,(jjd+1)/2,(jmd+1)/2,JJ,MMJ) ! J12,M12 are NOT doubled 

!             if(abs(cg1-cb1).gt.CHOP) stop 'Error: in CG_Save '

           DO IT=0,HO%tmax
            phasJT=(-1.d0)**(JJ+IT)
!           ........................ normalization factor
            AN_AB=sqrt(1.d0-delta_ab*phasJT)/(1.d0+delta_ab)
            AN_CD=sqrt(1.d0-delta_cd*phasJT)/(1.d0+delta_cd)
            if(AN_AB.le.1e-15) cycle
            if(AN_CD.le.1e-15) cycle
            AN_INV=1.d0/(AN_AB*AN_CD)
            DO MT=0,2*IT
             MMT=-IT+MT

!              call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)
!              call CJJ(1,1,2*IT,jtc,jtd,2*MMT,cb4)

              cb3  = CG_Save(1,(jta+1)/2,1,(jtb+1)/2,IT,MMT) ! J12,M12 are NOT doubled 
              cb4  = CG_Save(1,(jtc+1)/2,1,(jtd+1)/2,IT,MMT) ! J12,M12 are NOT doubled 

!             if(abs(cg1-cb3).gt.CHOP) stop 'Error: in CG_Save '

             sumando=AN_INV*cb1*cb2*cb3*cb4*                   &
     &        V_AB_CD_JTMT(jacoup,jbcoup,jccoup,jdcoup,JJ,IT,MMT)
             suma=suma+sumando
            END DO
           END DO

          END DO
          END DO
          
           H%ME2BM(iabcd) = suma
!           write(99) suma
       enddo ! iabcd         
!$OMP END DO 
!$OMP END PARALLEL
!     ..............................................................................
       if(lpr) write(*,*) '    Writing ME2B in m-scheme into file: '
       do iabcd = 1, H%iabcd_max
!         if(mod(iabcd,icheck).eq.0) call progress(iabcd/icheck)
         write(99) H%ME2BM(iabcd) ! suma
       enddo ! iabcd         
!
!     ..............................................................................
!       call cpu_time(finish)
!       print '("Time = ",f10.3," seconds.")',finish-start
!       if(allocated(CG_Save)) deallocate(CG_Save)
!     ..................................
!     read ME2B in m-scheme




!     ..................................
      else if(Input%IntJT.eq.0) then

      open(1,file=trim(INT_DIR)//trim(File%f2bm),form='unformatted',status='old')
      if(lpr) write(*,*) '    from ',trim(INT_DIR)//trim(File%f2bm)
      do iabcd=1,H%iabcd_max
           read(1) H%ME2BM(iabcd)
       enddo ! iabcd
      close(1)




!     ..................................
!     transform ME2B in J-scheme to JT-Scheme
!     ..................................
      elseif(Input%IntJT.eq.2) then
       if(lpr) &
       write(*,*) 'Main: read and transform H in J-scheme to JT-scheme'
       call GenerateV_JTMT()
!       if(Input%InME3B.eq.1) call V3B2V2B_JTMT()
       stop '  The ME2B in JT scheme is saved ...'





      elseif(Input%IntJT.eq.1) then
!     ..................................

!       open(99,file=trim(INT_DIR)//trim(File%f2bm),&
!               form='unformatted',status='unknown')

!        call  ME2J_Base()
         !open(99,file=trim(INT_DIR)//trim(File%f2bm),form='unformatted',status='unknown')
         !open(99,file='me2b.dat',form='unformatted',status='unknown')
!        call  ME2JT2M(jindex)

       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))

       open(48,file='fort.48',status='old')
11     read(48,*,end=12) ia,ib,ic,id,iJ,iT,MMT,ame_coup

       !print *, ia,ib,ic,id,iJ,iT,MMT,ame_coup
       hja=SPB%twoj(ia)/2.d0
       hjb=SPB%twoj(ib)/2.d0
       hjc=SPB%twoj(ic)/2.d0
       hjd=SPB%twoj(id)/2.d0
       ajtot=1.d0*iJ ! jtot(iJ)
       attot=1.d0*iT ! jttot(iT)

       phasab=(-1.d0)**(hja+hjb+ajtot+attot)
       phascd=(-1.d0)**(hjc+hjd+ajtot+attot)

       V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,ic,id,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,id,ic,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ia,ib,id,ic,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
!
       V_AB_CD_JTMT(ic,id,ia,ib,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ia,ib,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ib,ia,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ic,id,ib,ia,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme

       goto 11 
 12   continue
      close(48)

      !computing the two-body matrix elements in uncoupled basis

!     ..............................................................................
!      OMP
!     ..............................................................................
      icheck = H%iabcd_max/100

      print *, ' decoupling ...'
      print *, ' H%iabcd_max= ',H%iabcd_max
!$OMP PARALLEL DEFAULT(PRIVATE)
!SHARED(icheck,CG_Save,H,HO,tnljm,jindex,V_AB_CD_JTMT) 
!$OMP DO schedule(dynamic) 
      do iabcd = 1, H%iabcd_max
        ia = H%ka(iabcd)
        jja=tnljm%twoj(ia)  ! jang(ia)
        jma=tnljm%twom(ia)
        jta=tnljm%t(ia)      !mtisos(ia) 
        jna=tnljm%n(ia)
        jlja=tnljm%lj(ia)
        jacoup=jindex(ia)

!        print*,ia
        ib = H%kb(iabcd)
        jjb=tnljm%twoj(ib) ! jang(ib)
        jmb=tnljm%twom(ib)  ! mjang(ib)
        jtb=tnljm%t(ib)  !mtisos(ib)
        jnb=tnljm%n(ib)
        jljb=tnljm%lj(ib)
        jbcoup=jindex(ib)
        delta_ab=0.d0
        if(jna.eq.jnb .and. jlja.eq.jljb) delta_ab=1.d0
!        do ic=1,NLEV
        ic = H%kc(iabcd)
         jjc=tnljm%twoj(ic)   ! jang(ic)
         jmc=tnljm%twom(ic)   ! mjang(ic)
         jtc=tnljm%t(ic)      ! mtisos(ic)
         jnc=tnljm%n(ic)
         jljc=tnljm%lj(ic)
         jccoup=jindex(ic)
          id = H%kd(iabcd)
          jjd=tnljm%twoj(id)   !jang(id)
          jmd=tnljm%twom(id)   !mjang(id)
          jtd=tnljm%t(id)      !mtisos(id)
          jnd=tnljm%n(id)
          jljd=tnljm%lj(id)
          jdcoup=jindex(id)
          delta_cd=0.d0
        if(jnc.eq.jnd .and. jljc.eq.jljd) delta_cd=1.d0

          suma=0.d0
          DO JJ=0,HO%twojmax
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ

!              call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)
!              call CJJ(jjc,jjd,2*JJ,jmc,jmd,2*MMJ,cb2)
              cb1  = CG_Save((jja+1)/2,(jma+1)/2,(jjb+1)/2,(jmb+1)/2,JJ,MMJ) 
              cb2  = CG_Save((jjc+1)/2,(jmc+1)/2,(jjd+1)/2,(jmd+1)/2,JJ,MMJ) 

!             if(abs(cg1-cb1).gt.CHOP) stop 'Error: in CG_Save '

           DO IT=0,HO%tmax
            phasJT=(-1.d0)**(JJ+IT)
!           ........................ normalization factor
            AN_AB=sqrt(1.d0-delta_ab*phasJT)/(1.d0+delta_ab)
            AN_CD=sqrt(1.d0-delta_cd*phasJT)/(1.d0+delta_cd)
            if(AN_AB.le.1e-15) cycle
            if(AN_CD.le.1e-15) cycle
            AN_INV=1.d0/(AN_AB*AN_CD)
            DO MT=0,2*IT
             MMT=-IT+MT

!              call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)
!              call CJJ(1,1,2*IT,jtc,jtd,2*MMT,cb4)

              cb3  = CG_Save(1,(jta+1)/2,1,(jtb+1)/2,IT,MMT) ! J12,M12
              cb4  = CG_Save(1,(jtc+1)/2,1,(jtd+1)/2,IT,MMT) ! J12,M12

!             if(abs(cg1-cb3).gt.CHOP) stop 'Error: in CG_Save '

             sumando=AN_INV*cb1*cb2*cb3*cb4*                   &
     &        V_AB_CD_JTMT(jacoup,jbcoup,jccoup,jdcoup,JJ,IT,MMT)
             suma=suma+sumando
            END DO
           END DO

          END DO
          END DO
          
           H%ME2BM(iabcd) = suma
       enddo ! iabcd         
!$OMP END DO 
!$OMP END PARALLEL
!     ..............................................................................
          if(lpr) write(*,*) '    Writing ME2B in m-scheme into file: '
          do iabcd = 1, H%iabcd_max
!           if(mod(iabcd,icheck).eq.0) call progress(iabcd/icheck)
             write(911) H%ME2BM(iabcd) ! suma
          enddo ! iabcd         
!
!     ................................................
!     Read ME2B in J and transform them into M-scheme
!     ...............................................
      elseif(Input%IntJT.eq.4) then
             open(99,file=trim(INT_DIR)//trim(File%f2bm),&
                     form='unformatted',status='unknown')
             call GenerateV_J2M()
!     ..................................
         
      else
            if(lpr) print*,'change interaction option'
            stop
      end if
!  ...................
      
      return

      END


