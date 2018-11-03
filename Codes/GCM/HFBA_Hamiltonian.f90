        subroutine Hamiltonian() 
!     ..........................
        USE VAPHFB_PAR
        implicit none
        integer mit,nn
        real*8, dimension(:,:,:,:,:,:,:), allocatable :: V_AB_CD_JTMT
        real*8, dimension(:,:,:,:,:,:), allocatable :: V_AB_CD_J,V_AB_CD_JT

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


      write(*,*) ' -----------------------'
      write(*,*) ' idx   n+1    l    2j '
      do nlj=1,HO%nljmax
!         nlshort(nlj)=SPB%n(nlj)+1   ! 1d5/2; n =1,2,..
!         lshort(nlj)=SPB%l(nlj)
!         jshort(nlj)=SPB%twoj(nlj)
         write(*,'(4i5)') nlj,SPB%n(nlj)+1,SPB%l(nlj),SPB%twoj(nlj)
      enddo    

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
      if(kk_n .ne. HO%NLEV) stop 'Main: NLEV is wrong ...'
!     ............... quantum numbers for each orbital
      if(.NOT. ALLOCATED(tnljm%t)) ALLOCATE(tnljm%t(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%n)) ALLOCATE(tnljm%n(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%l)) ALLOCATE(tnljm%l(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%lj)) ALLOCATE(tnljm%lj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twoj)) ALLOCATE(tnljm%twoj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twom)) ALLOCATE(tnljm%twom(1:HO%NLEV))
      if(.NOT. ALLOCATED(HO%tts)) ALLOCATE(HO%tts(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%level)) ALLOCATE(tnljm%level(-HO%twojmax:HO%twojmax,0:HO%ljmax,0:HO%nmax,0:1))
      write(*,*) '.................................'
      write(*,*) '| Single-Particle Basis         |' 
      write(*,*) '.................................'
      write(*,*)  'k   t   n+1   l   twoj   twom'
      do kk=1,kk_n       ! max number of s.p. state
          tnljm%t(kk) = mtisos(kk) 
          tnljm%n(kk) = nlindex(kk)
          tnljm%l(kk) = lang(kk)
          tnljm%twoj(kk) = jang(kk)
          tnljm%lj(kk)   = (2*lang(kk)+jang(kk)-1)/2
          tnljm%twom(kk) = mjang(kk)
          tnljm%level(mjang(kk),tnljm%lj(kk),tnljm%n(kk),mit(mtisos(kk))) = kk
!          write(*,'(5i4)')  kk,mjang(kk),tnljm%lj(kk),tnljm%n(kk),mit(mtisos(kk))
          write(*,20)  kk,tnljm%t(kk),tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),'/2',tnljm%twom(kk),'/2'

! ............
          nn = 2*tnljm%n(kk)+tnljm%l(kk)
          if (nn.lt.10) then
              write(HO%tts(kk),19) tnljm%n(kk),HO%tl(tnljm%l(kk)), &
     &                 tnljm%twoj(kk),'/2',tnljm%twom(kk),'/2',    &
     &                 HO%tp(mip(iv(tnljm%l(kk)))) 
               else
                  write(HO%tts(kk),21) tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),tnljm%twom(kk)
               endif

!          write(*,'(i3,a14)') kk, HO%tts(kk)
      end do
19   format(1h[,i2,a,i2,a,i2,a,1h],a)
20   format(5i4,a,i4,a)
21   format(4i2)
!   .................................
!      one-body energy
       if(.NOT. ALLOCATED(H%ME1BM)) ALLOCATE(H%ME1BM(1:HO%NLEV,1:HO%NLEV))
       CALL kinetic(H%ME1BM,H%ddd_kin)

!     determine the dimension
      write(*,*) '........... prepare the two-body matrix elements ...'
      iabcd = 0
      do ia=1,HO%NLEV
      do ib=ia,HO%NLEV
         do ic=1,HO%NLEV
         do id=1,HO%NLEV
         if(tnljm%t(ia)+tnljm%t(ib) .ne. tnljm%t(ic)+tnljm%t(id))             cycle   ! isospin conservation
         if(tnljm%twom(ia)+tnljm%twom(ib) .ne. tnljm%twom(ic)+tnljm%twom(id)) cycle   ! m conservation
         if(iv(tnljm%l(ia)+tnljm%l(ib)) .ne. iv(tnljm%l(ic)+tnljm%l(id)))     cycle   ! parity
!         if(ia+ib .gt. ic+id)                                                 cycle
         iabcd     = iabcd + 1
      enddo
      enddo
      enddo
      enddo

      H%iabcd_max = iabcd
      write(*,*) 'No. of stored two-body matrix elements=',H%iabcd_max

      if(.NOT. ALLOCATED(H%ka))  ALLOCATE(H%ka(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kb))  ALLOCATE(H%kb(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kc))  ALLOCATE(H%kc(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%kd))  ALLOCATE(H%kd(1:H%iabcd_max))
      if(.NOT. ALLOCATED(H%ME2BM))  ALLOCATE(H%ME2BM(1:H%iabcd_max))
       write(*,'(a,f6.3,a)')'Main: allocated ka,kb,kc,kd with memory:',&
     & 4*sizeof(H%ka)/(1024*1024*1024.0),'G byte ..'
      write(*,'(a,f6.3,a)')'Main: allocated H%ME2B with memory:', &
     & sizeof(H%ME2BM)/(1024*1024*1024.0),'G byte ..'


      iabcd = 0
      do ia=1,HO%NLEV
      do ib=ia,HO%NLEV
         do ic=1,HO%NLEV
         do id=1,HO%NLEV
         if(tnljm%t(ia)+tnljm%t(ib) .ne. tnljm%t(ic)+tnljm%t(id))             cycle   ! isospin conservation
         if(tnljm%twom(ia)+tnljm%twom(ib) .ne. tnljm%twom(ic)+tnljm%twom(id)) cycle   ! m conservation
         if(iv(tnljm%l(ia)+tnljm%l(ib)) .ne. iv(tnljm%l(ic)+tnljm%l(id)))     cycle   ! parity
!         if(ia+ib .gt. ic+id)                                                 cycle
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


      !Initializing the two-body matrix elements in coupled basis
       if(.NOT. ALLOCATED(V_AB_CD_J))  &
     & ALLOCATE(V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5))
       if(.NOT. ALLOCATED(V_AB_CD_JT))  &
     & ALLOCATE(V_AB_CD_JT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1))
      do ia=1,HO%nljmax
       do ib=1,HO%nljmax
        do ic=1,HO%nljmax
         do id=1,HO%nljmax
          do iJ=0,HO%twojmax
!  ...................
           do iT=0,HO%tmax
              V_AB_CD_JT(ia,ib,ic,id,iJ,iT) =0.d0
           end do
           do it12=0,5
            V_AB_CD_J(ia,ib,ic,id,iJ,it12)=0.d0
           end do
!  ...................
          end do
         end do
        end do
       end do
      end do

       open(99,file='../Int/IMSRG_'//Input%cFlow//'_'//Input%cIntID &
     &  //'_m1m2m3m4_COM1B.dat',form='unformatted',status='unknown')
!   ......reading the two-body matrix elements in coupled basis from file
      write(*,*) 'Main: read hamiltonian matrix elements'
      if(Input%IntJT.eq.1) then
!  .........................................................
!       if(Input%IntIMSRG .ne. 1) then
!         open(1,file='fort.88',status='old')
!         write(*,*) ' from ',Input%cIntID,'.dat'
!       endif
!     ....................................................
!      READ H with Isospin Symmetry
!     ..................................................
       if(Input%IntIMSRG .eq. 1 .and. Input%IsoSpin .ne.0) then
        open(1,file='../Int/IMSRG_'//Input%cFlow//'_'//Input%cIntID &
     &                //'_JT.dat',status='old')
        write(*,*) ' from IMSRG_',Input%cFlow,'_',Input%cIntID,'_JT.dat'
4      read(1,*,end=5) ia,ib,ic,id,iJ,iT,ame_coup
       hja=SPB%twoj(ia)/2.d0
       hjb=SPB%twoj(ib)/2.d0
       hjc=SPB%twoj(ic)/2.d0
       hjd=SPB%twoj(id)/2.d0

       ajtot=1.d0*iJ ! jtot(iJ)
       attot=1.d0*iT ! jttot(iT)
       phasab=(-1.d0)**(hja+hjb+ajtot+attot)
       phascd=(-1.d0)**(hjc+hjd+ajtot+attot)

       V_AB_CD_JT(ia,ib,ic,id,iJ,iT)= ame_coup*H%ddd_tbme
       V_AB_CD_JT(ib,ia,ic,id,iJ,iT)= ame_coup*phasab*H%ddd_tbme
       V_AB_CD_JT(ib,ia,id,ic,iJ,iT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JT(ia,ib,id,ic,iJ,iT)= ame_coup*phascd*H%ddd_tbme
!
       V_AB_CD_JT(ic,id,ia,ib,iJ,iT)= ame_coup*H%ddd_tbme
       V_AB_CD_JT(id,ic,ia,ib,iJ,iT)= ame_coup*phascd*H%ddd_tbme
       V_AB_CD_JT(id,ic,ib,ia,iJ,iT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JT(ic,id,ib,ia,iJ,iT)= ame_coup*phasab*H%ddd_tbme

       goto 4
  5   continue
      close(1)

       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))
      do iT=0,HO%tmax
      do MMT=0,HO%tmax
      do ia=1,HO%nljmax
       do ib=1,HO%nljmax
        do ic=1,HO%nljmax
         do id=1,HO%nljmax
          do iJ=0,HO%twojmax
!  ...................
            V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT) = V_AB_CD_JT(ia,ib,ic,id,iJ,iT)
          end do
          end do
         end do
        end do
       end do
      end do
      end do
!  ...................
      endif
!  ...................................
!      READ H without Isospin Symmetry
!  ....................................
      if(Input%IntIMSRG .eq. 1 .and. Input%IsoSpin .eq.0) then

       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))
      write(*,*) ' from IMSRG_',Input%cFlow,'_',Input%cIntID,'_JTMT.dat'
       open(1,file='../Int/IMSRG_'//Input%cFlow//'_'//Input%cIntID//&
     &  '_JTMT.dat',status='old')
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
      endif

      !computing the two-body matrix elements in uncoupled basis
      write(*,*) 'Main: uncouple hamiltonian matrix elements'

      icheck = H%iabcd_max/10

      do iabcd = 1, H%iabcd_max

!       if(mod(iabcd,icheck) .eq. 0) write(*,*) '........ complete',iabcd/icheck,'%.......'
        if(mod(iabcd,icheck).eq.0) call progress(iabcd/icheck)
!      do ia=1,NLEV
        ia = H%ka(iabcd)
        jja=tnljm%twoj(ia)  ! jang(ia)
        jma=tnljm%twom(ia)
        jta=tnljm%t(ia)      !mtisos(ia) 
        jna=tnljm%n(ia)
        jlja=tnljm%lj(ia)
        jacoup=jindex(ia)

!        print*,ia
!       do ib=1,NLEV
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
         jjc=jang(ic)
         jmc=mjang(ic)
         jtc=mtisos(ic)
         jnc=tnljm%n(ic)
         jljc=tnljm%lj(ic)
         jccoup=jindex(ic)
!         do id=1,NLEV
          id = H%kd(iabcd)
          jjd=jang(id)
          jmd=mjang(id)
          jtd=mtisos(id)
          jnd=tnljm%n(id)
          jljd=tnljm%lj(id)
          jdcoup=jindex(id)
          delta_cd=0.d0
        if(jnc.eq.jnd .and. jljc.eq.jljd) delta_cd=1.d0

          suma=0.d0
          DO JJ=0,HO%twojmax
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ
              call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)
              call CJJ(jjc,jjd,2*JJ,jmc,jmd,2*MMJ,cb2)
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
             call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)
             call CJJ(1,1,2*IT,jtc,jtd,2*MMT,cb4)
             sumando=AN_INV*cb1*cb2*cb3*cb4*                   &
     &        V_AB_CD_JTMT(jacoup,jbcoup,jccoup,jdcoup,JJ,IT,MMT)
             suma=suma+sumando
            END DO
           END DO

          END DO
          END DO
          
           H%ME2BM(iabcd) = suma
           write(99) suma
!         end do
!        end do
!       end do
!      end do
       enddo ! iabcd         

      else if(Input%IntJT.eq.0) then

      open(1,file='../Int/IMSRG_'//Input%cFlow//'_'//Input%cIntID//  &
     & '_m1m2m3m4_COM1B.dat',form='unformatted',status='old')
      write(*,*) ' from IMSRG_',Input%cFlow,'_', &
     & Input%cIntID,'_m1m2m3m4_COM1B.dat'

      do iabcd=1,H%iabcd_max
           read(1) H%ME2BM(iabcd)
       enddo ! iabcd

      elseif(Input%IntJT.eq.2) then
       write(*,*) 'Main: read and transform H in J-scheme to JT-scheme'

       call GenerateV_JTMT

      else
       print*,'change interaction option'
       stop
      end if
!  ...................
      
      return

      END


