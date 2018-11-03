        subroutine Hamiltonian(lpr) 
!     ..........................
        USE VAPHFB_PAR
        use omp_lib
        implicit none
        logical lpr
        integer mitf,nn
        real*8, dimension(:,:,:,:,:,:,:), allocatable :: V_AB_CD_JTMT
        !real*8, dimension(:,:,:,:,:,:), allocatable :: V_AB_CD_J,V_AB_CD_JT
        real*8 start,finish
        integer ii,kk,mj,kk_n
        integer ia,ib,ic,id,na,lja,nb,ljb,nc,ljc,nd,ljd
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

!      ..... one-body
!          File%f1b="../Int/IMSRG_"//trim(Input%cFlow)//"_SPE_"//trim(Input%cIntID)//&
!          & "_"//trim(Input%cValID)//'_'//Input%chwHO//".dat"
!

!      ..... two-body

           File%f2bJTMT='IMSRG_'//trim(Input%cFlow)//'_'//trim(Input%cIntID)//&
        &  Input%ctpp//trim(Input%cValID)//'_'//Input%chwHO//'_JTMT.dat'

           File%f2bm='IMSRG_'//trim(Input%cFlow)//'_'//trim(Input%cIntID)//&
         &Input%ctpp//trim(Input%cValID)//'_'//Input%chwHO//&
         &'_m1m2m3m4.dat'

        print *,' me2b(JTMT) file:',File%f2bJTMT 
        print *,' me2b(m) file:',File%f2bm

!      ....... remove spaces in the filenames
      
       call StripSpaces(File%f1b)
       call StripSpaces(File%f2bm)
       call StripSpaces(File%f2bJTMT)

       if(Input%IntJT.eq.0) INT_DIR = find_file("GCM_ME_FILES",trim(File%f2bm))
      
!       if(Input%IntJT .ne.0) then
          open(41,file='chi2b3b.int',status='old')
          read(41,'(a14,a)') ctemp,File%IMSRG_Hme1b
          read(41,'(a14,a)') ctemp,File%IMSRG_Hme2b   ! in J-scheme
          close(41)
          INT_DIR = find_file("GCM_ME_FILES",File%IMSRG_Hme1b)
        File%f2bJ = File%IMSRG_Hme2b
        print *,' me2b(J) file:',File%IMSRG_Hme2b 
        print *,' me1b file:',File%IMSRG_Hme1b 
          
        print *,'  matrix elements files are found in: ',trim(INT_DIR)
!       endif

       print *, ' ======================== '
       print *, ' Details of Calculations  '
       print *, ' ======================== '

      if(lpr) then
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
      if(.NOT. ALLOCATED(tnljm%nlj)) ALLOCATE(tnljm%nlj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twoj)) ALLOCATE(tnljm%twoj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twom)) ALLOCATE(tnljm%twom(1:HO%NLEV))
      if(.NOT. ALLOCATED(HO%tts)) ALLOCATE(HO%tts(1:HO%NLEV))
       if(.NOT. ALLOCATED(tnljm%level)) &
     & ALLOCATE(tnljm%level(-HO%twojmax:HO%twojmax,0:HO%ljmax,0:HO%nmax,0:1))

      if(lpr) then
        write(*,*) '.................................'
        write(*,*) '| Single-Particle Basis         |' 
        write(*,*) '.................................'
        write(*,*)  'k   t   n+1   l   twoj   twom'
      endif
      do kk=1,kk_n       ! max number of s.p. state
          tnljm%t(kk) = mtisos(kk) 
          tnljm%n(kk) = nlindex(kk)
          tnljm%l(kk) = lang(kk)
          tnljm%twoj(kk) = jang(kk)
          tnljm%lj(kk)   = (2*lang(kk)+jang(kk)-1)/2
          tnljm%nlj(kk) = jindex(kk)
          tnljm%twom(kk) = mjang(kk)
         tnljm%level(mjang(kk),tnljm%lj(kk),tnljm%n(kk),mitf(mtisos(kk))) = kk
     !  if(lpr) &   
     ! &write(*,20)  kk,tnljm%t(kk),tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),'/2',tnljm%twom(kk),'/2'

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
       print *, ' ...............'
       print *, ' -> Hamiltonian '

!      one-body energy
       if(.NOT. ALLOCATED(H%ME1BM)) ALLOCATE(H%ME1BM(1:HO%NLEV,1:HO%NLEV))

!     ............... read ME1B from files 
      !if(Input%IntJT.eq.3) call Generate_ME1B(H%E0,H%ME1BM)   ! read ME1B from the files by Nathan
      !if(Input%IntJT.ne.3) 

      call kinetic(H%ME1BM,H%ddd_kin)
 
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
!      if(lpr) then
       write(*,*) '   No. of stored two-body matrix elements:',H%iabcd_max
       write(*,'(a,f6.3,a)')'    Main: allocated ka,kb,kc,kd with &
     & memory:', 4*sizeof(H%ka)/(1024*1024*1024.0),'G byte ..'
       write(*,'(a,f6.3,a)')'    Main: allocated H%ME2B with memory:', &
     & sizeof(H%ME2BM)/(1024*1024*1024.0),'G byte ..'
!      endif

!      if(lpr)
      print *, ' -> Initializing ka,kb,kc,kd ...'
      iabcd = 0
      icheck = HO%NLEV/100
      if(icheck.eq.0) icheck=HO%NLEV/10
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(iv,HO,tnljm,H) REDUCTION(+: iabcd)
!!$OMP DO SCHEDULE(DYNAMIC) 
      do ia=1,HO%NLEV
!        if(mod(ia,icheck).eq.0) call progress(ia/icheck)
        if(mod(ia,icheck).eq.0) call progress_bar(int(ia/icheck))  
      do ib=ia+1,HO%NLEV
         do ic=1,HO%NLEV
         do id=ic+1,HO%NLEV
         if(tnljm%t(ia)+tnljm%t(ib) .ne. tnljm%t(ic)+tnljm%t(id))             cycle   ! isospin conservation
         if(tnljm%twom(ia)+tnljm%twom(ib) .ne. tnljm%twom(ic)+tnljm%twom(id)) cycle   ! m conservation
         if(iv(tnljm%l(ia)+tnljm%l(ib)) .ne. iv(tnljm%l(ic)+tnljm%l(id)))     cycle   ! parity
         if(ia+ib .gt. ic+id)                                                 cycle
         iabcd     = iabcd + 1
         H%ka(iabcd) = ia
         H%kb(iabcd) = ib
         H%kc(iabcd) = ic
         H%kd(iabcd) = id
         H%ME2BM(iabcd) = 0.d0
      enddo
      enddo
      enddo
      enddo
!!$OMP END DO
!!$OMP END PARALLEL
      write(*,*)

!     ..................................................
!     decoupling ME2B from JTMT to m-scheme
!     ..................................................

      if(Input%IntJT.eq.0) then
!     ..............................................................................
!     read the two-body matrix elements in m-scheme directly
!     ..............................................................................
          open(1,file=trim(INT_DIR)//trim(File%f2bm),form='unformatted',status='old')
          write(*,*) '    from ',trim(INT_DIR)//trim(File%f2bm)
          do iabcd=1,H%iabcd_max
             read(1) H%ME2BM(iabcd)
          enddo ! iabcd
          return


      else if(Input%IntJT.eq.1) then
!     ..............................................................................
!     read the two-body matrix elements in J-scheme and transform into m-scheme 
!     ..............................................................................
          call INT_ME2B_J2M()
          return


       else if(Input%IntJT.eq.2) then
!     ..............................................................................
!      read the two-body matrix elements in J-scheme and transform into
!      JT and finally to m-scheme
!     ..............................................................................
           write(*,*) 'Main: transform H in J-scheme to JT-scheme'
           call INT_ME2B_J2JT()
           !if(Input%InME3B.eq.1) call V3B2V2B_JTMT()
           write(*,*) 'Main: transform H in JT-scheme to m-scheme'
           call INT_ME2B_JT2M()
           return


      else if(Input%IntJT.eq.3) then
!     ..............................................................................
!     transform ME2B in J-scheme to JT-Scheme
!     ..................................
          write(*,*) 'Main: extract H in JT-scheme from files directly'

          File%NN   ='chi2b_srg0953_eMax06_hwHO020.me2j'
          File%p1p2 ='tpp_eMax06.me2j'
          call  INT_ME2J_JT()     

          write(*,*) 'Main: transform H in JT-scheme to m-scheme'
          call INT_ME2B_JT2M()
          return
!     ................................................
!       call EOM_ME2B_J2M()
!       write(*,*) 'Main: read and transform H (by Nathan) &
!      & in J-scheme to JT-scheme'
!       call Generate_ME2B()
!     ................................................
!     Read ME2B in J and transform them into M-scheme
!     ...............................................
         
      else
       print*,'change interaction option'
       stop
      end if
!  ...................
      
      return

      END


