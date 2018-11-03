        subroutine TBME4DBD(lpr) 
!     ..........................
        USE VAPHFB_PAR
        implicit none

        logical lpr
        integer mitf,ii,kk,mj,kk_n
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
        integer icheck,k1,k2,k3,k4
        real*8  DBD_GT,DBD_FM,DBD_TE 
        integer name_emax1,name_emax2,hw1,hw2
        integer n1,n2,n3,n4,nlj1,lj1,lj2,lj3,lj4,nlj2,nlj3,nlj4

        Integer pf_lj(1:4),pf_n(1:4)
        character(500) find_file


      if(lpr) write(*,*) ' -----------------------'
      if(lpr) write(*,*) ' idx   n+1    l    2j '
      do nlj=1,HO%nljmax
         if(lpr) write(*,'(4i5)') nlj,SPB%n(nlj)+1,SPB%l(nlj),SPB%twoj(nlj)
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
      if(kk_n .ne. HO%NLEV) stop 'DBD_TBME: NLEV is wrong ...'
!     ............... quantum numbers for each orbital
      if(.NOT. ALLOCATED(tnljm%t)) ALLOCATE(tnljm%t(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%n)) ALLOCATE(tnljm%n(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%l)) ALLOCATE(tnljm%l(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%lj)) ALLOCATE(tnljm%lj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twoj)) ALLOCATE(tnljm%twoj(1:HO%NLEV))
      if(.NOT. ALLOCATED(tnljm%twom)) ALLOCATE(tnljm%twom(1:HO%NLEV))
       if(.NOT. ALLOCATED(tnljm%level)) &
     & ALLOCATE(tnljm%level(-HO%twojmax : HO%twojmax,0:HO%ljmax,0:HO%nmax,0:1))

      if(lpr) then
        write(*,*) '.................................'
        write(*,*) '| Single-Particle Basis         |' 
        write(*,*) '.................................'
        write(*,*)  '  k  t  n+1   l   twoj   twom'
      endif
       do kk=1,kk_n       ! max number of s.p. state
          tnljm%t(kk) = mtisos(kk) 
          tnljm%n(kk) = nlindex(kk) ! starting from 1
          tnljm%l(kk) = lang(kk)
          tnljm%twoj(kk) = jang(kk)
          tnljm%lj(kk)   = (2*lang(kk)+jang(kk)-1)/2
          tnljm%twom(kk) = mjang(kk)
          tnljm%level(mjang(kk),tnljm%lj(kk),tnljm%n(kk),mitf(mtisos(kk))) = kk

        if(lpr) write(*,20)  kk,tnljm%t(kk),tnljm%n(kk),tnljm%l(kk),tnljm%twoj(kk),'/2',tnljm%twom(kk),'/2'
      end do
20   format(5i4,a,i4,a)

!   .................................
! ....................................
!  TWO-BODY Transition Matrix Elements 
! ....................................
!
! ................... initialization for the DBD operators

      if(.NOT. ALLOCATED(DBD%FM2B)) ALLOCATE(DBD%FM2B(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%Twojmax))
      if(.NOT. ALLOCATED(DBD%GT2B)) ALLOCATE(DBD%GT2B(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%Twojmax))
      if(.NOT. ALLOCATED(DBD%TE2B)) ALLOCATE(DBD%TE2B(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%Twojmax))

      write(*,*) '  ->Main: reading DBD 2BME in J-scheme'

      print *, ' Input%NME_Type=',Input%NME_Type
!     ............................... read from Jon
      if(Input%NME_Type.eq.0 .and. Input%cIntID(1:3) == 'KB3') then
!     ............. model space
!     1: 0f7/2
!     2: 1p3/2
!     3: 0f5/2
!     4: 1p1/2
!     ....................
       pf_lj(1) = 6; pf_lj(2) = 2; pf_lj(3) = 5; pf_lj(4)=1
       pf_n(1)  = 0; pf_n(2)  = 1; pf_n(3)  = 0; pf_n(4) =1

        if(trim(Input%cFlow) =='s000') then
         File%GT='Jon_0nuME_48Ca_pf_hwHO011_GT.dat'
         File%FM='Jon_0nuME_48Ca_pf_hwHO011_FM.dat'
        else
         File%GT='Jon_0nuME_48Ca_pf_'//trim(Input%cIntID)//'hwHO011_GT.dat'
         File%FM='Jon_0nuME_48Ca_pf_'//trim(Input%cIntID)//'hwHO011_FM.dat'
        endif
      INT_DIR = find_file("DBD_ME_FILES",File%GT)

      open(10,file=trim(INT_DIR)//File%GT,status='old')
      write(*,*) ' GT part from ',trim(INT_DIR)//File%GT
 31   read(10,*,end=41) k1,k2,k3,k4,JJ,DBD_GT

        if(JJ.gt. HO%TwoJmax) stop 'Main: TwoJMax is too small'
!     ................
!     only for pf
!     ................
      n1 = pf_n(k1) 
      n2 = pf_n(k2) 
      n3 = pf_n(k3) 
      n4 = pf_n(k4) 
       lj1 = pf_lj(k1)
       lj2 = pf_lj(k2)
       lj3 = pf_lj(k3)
       lj4 = pf_lj(k4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)     ! n1 starts from 0 
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%GT2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_GT
       DBD%GT2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_GT
       write(20,'(4i5,i6,f12.8)') nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT
      goto 31
 41   continue
      close(10)

      write(*,*) ' FM part from ',trim(INT_DIR)//File%FM
      open(10,file=trim(INT_DIR)//File%FM,status='old')
131   read(10,*,end=141) k1,k2,k3,k4,JJ,DBD_FM
      if(JJ.gt. HO%TwoJmax) stop 'Main: TwoJMax is too small'
!     ................
!     only for pf
!     ................
      n1 = pf_n(k1)
      n2 = pf_n(k2)
      n3 = pf_n(k3)
      n4 = pf_n(k4)
       lj1 = pf_lj(k1)
       lj2 = pf_lj(k2)
       lj3 = pf_lj(k3)
       lj4 = pf_lj(k4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)     ! n1 starts from 0 
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%FM2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_FM
       DBD%FM2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_FM
      goto 131
 141   continue
      close(10)

      else if(Input%NME_Type.eq.4) then
         write(*,*) '   from the file', &
     &             ' ../../NME/Mihai_0nuME_76Ge_pf5g9_hwHO09.dat.'

      open(10,file="../../NME/Mihai_0nuME_76Ge_pf5g9_hwHO09.dat",&
     &              status='old')
 231   read(10,*,end=241) nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT,DBD_FM,DBD_TE
       pf_lj(1) = 5; pf_lj(2) = 2; pf_lj(3) = 1; pf_lj(4)=8
       pf_n(1)  = 0; pf_n(2)  = 1; pf_n(3)  = 1; pf_n(4) =0

       n1 = pf_n(nlj1)
       n2 = pf_n(nlj2)
       n3 = pf_n(nlj3)
       n4 = pf_n(nlj4)
       lj1 = pf_lj(nlj1)
       lj2 = pf_lj(nlj2)
       lj3 = pf_lj(nlj3)
       lj4 = pf_lj(nlj4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%GT2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_GT
       DBD%GT2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_GT
       DBD%FM2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_FM
       DBD%FM2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_FM
       DBD%TE2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_TE
       DBD%TE2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_TE
       write(30,'(4i5,i6,f12.8)') nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT
      goto 231
  241   continue
      close(10) 

      else if(Input%NME_Type.eq.1 .and. Input%cIntID == 'KB3') then
      write(*,*) '  from the file ../NME/Jon_0nuME_48Ca_pf_hwHO011_GT.dat ....'
!    ............................... read from Jon
       pf_lj(1) = 6; pf_lj(2) = 2; pf_lj(3) = 5; pf_lj(4)=1
       pf_n(1)  = 0; pf_n(2)  = 1; pf_n(3)  = 0; pf_n(4) =1
 
      open(10,file="../NME/Jon_0nuME_48Ca_pf_hwHO011_GT.dat",status='old')
 51   read(10,*,end=61) nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT
       n1 = pf_n(nlj1)
       n2 = pf_n(nlj2)
       n3 = pf_n(nlj3)
       n4 = pf_n(nlj4)
       lj1 = pf_lj(nlj1)
       lj2 = pf_lj(nlj2)
       lj3 = pf_lj(nlj3)
       lj4 = pf_lj(nlj4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%GT2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_GT
       DBD%GT2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_GT
       write(30,'(4i5,i6,f12.8)') nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT
      goto 51
 61   continue
      close(10)
      write(*,*) '  from the file ../NME/Jon_0nuME_48Ca_pf_hwHO011_FM.dat ....'
      open(11,file="../NME/Jon_0nuME_48Ca_pf_hwHO011_FM.dat",status='old')
151   read(11,*,end=161) nlj1,nlj2,nlj3,nlj4,JJ,DBD_FM
       n1 = pf_n(nlj1)
       n2 = pf_n(nlj2)
       n3 = pf_n(nlj3)
       n4 = pf_n(nlj4)
       lj1 = pf_lj(nlj1)
       lj2 = pf_lj(nlj2)
       lj3 = pf_lj(nlj3)
       lj4 = pf_lj(nlj4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%FM2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_FM
       DBD%FM2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_FM
      goto 151
161   continue
      close(11)

!   ...........................................................
      else if(Input%NME_Type.eq.2 .and. Input%cIntID == 'KB3') then
      write(*,*) '  from the file &
      &../NME/Mihai_0nuME_48Ca_pf_hwHO011.dat ....'
       pf_lj(1) = 6; pf_lj(2) = 2; pf_lj(3) = 5; pf_lj(4)=1
       pf_n(1)  = 0; pf_n(2)  = 1; pf_n(3)  = 0; pf_n(4) =1
      open(10,file="../NME/Mihai_0nuME_48Ca_pf_hwHO011.dat",status='old')
 71   read(10,*,end=81) nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT,DBD_FM,DBD_TE
       n1 = pf_n(nlj1)
       n2 = pf_n(nlj2)
       n3 = pf_n(nlj3)
       n4 = pf_n(nlj4)
       lj1 = pf_lj(nlj1)
       lj2 = pf_lj(nlj2)
       lj3 = pf_lj(nlj3)
       lj4 = pf_lj(nlj4)
!     ................
       nlj1 = SPB%nlj(n1,lj1)
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%GT2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_GT
       DBD%GT2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_GT
       DBD%FM2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_FM
       DBD%FM2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_FM
       DBD%TE2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_TE
       DBD%TE2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_TE
       write(30,'(4i5,i6,3f12.8)') nlj1,nlj2,nlj3,nlj4,JJ,DBD_GT,DBD_FM,DBD_TE
      goto 71
 81   continue
      close(10)

      else if(Input%NME_Type.eq.3) then

      name_emax1 =mod(HO%emax/10,10) + 48
      name_emax2 =mod(HO%emax,10) + 48
      hw1= mod(Input%ihwHO/10,10)+48
      hw2= mod(Input%ihwHO,10)+48
      if(trim(Input%cFlow) =='s000') then
         File%GT="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax0.0_GT.dat"
         File%FM="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax0.0_FM.dat"
         File%TE="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax0.0_TE.dat"

      else if(trim(Input%cFlow) =='s999') then
         File%GT="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax1.0_"//trim(Input%cIntID)              &
     &          //"_GT.dat"
         File%FM="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax1.0_"//trim(Input%cIntID)              &
     &          //"_FM.dat"
         File%TE="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax1.0_"//trim(Input%cIntID)              &
     &          //"_TE.dat"
       else
         File%GT="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax"//trim(Input%cFlow(2:4))                 &
     &           //"_"//trim(Input%cIntID)                   &
     &          //"_GT.dat"
         print *,'GT file:',File%GT
         File%FM="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax"//trim(Input%cFlow(2:4))                 &
     &           //"_"//trim(Input%cIntID)                   &
     &          //"_FM.dat"
         File%TE="DBD0nu"                                     &
     &          //'_eMax'//char(name_emax1)//char(name_emax2) &
     &          //'_lMax'//char(name_emax1)//char(name_emax2) &
     &          //'_hwHO0'//char(hw1)//char(hw2)               &
     &          //"_sMax"//trim(Input%cFlow(2:4))                 &
     &           //"_"//trim(Input%cIntID)                   &
     &          //"_TE.dat"

      endif

      INT_DIR = find_file("DBD_ME_FILES",File%GT)
      open(10,file=trim(INT_DIR)//File%GT,status='old')
      write(*,*) ' GT part from ',trim(INT_DIR)//File%GT
 91   read(10,*,end=101) n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,DBD_GT
!     ................
       nlj1 = SPB%nlj(n1,lj1)  ! SPB%nlj(n,lj), where n starts from 0
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%GT2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_GT
       DBD%GT2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_GT
      goto 91
101   continue
      close(10)

      open(110,file=trim(INT_DIR)//File%FM,status='old')
      write(*,*) ' Fermi part from ',trim(INT_DIR)//File%FM
191   read(110,*,end=201) n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,DBD_FM
!     ................
       nlj1 = SPB%nlj(n1,lj1)  ! SPB%nlj(n,lj), where n starts from 0
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%FM2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_FM
       DBD%FM2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_FM
!       write(300,'(4i3,5i4,f10.6)') n1-1,n2-1,n3-1,n4-1,lj1,lj2,lj3,lj4,JJ,DBD_FM
      goto 191
201   continue
      close(110)

      open(210,file=trim(INT_DIR)//File%TE,status='old')
      write(*,*) ' Tensor part from ',trim(INT_DIR)//File%TE
291   read(210,*,end=301) n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,DBD_TE
!     ................
       nlj1 = SPB%nlj(n1,lj1)  ! SPB%nlj(n,lj), where n starts from 0
       nlj2 = SPB%nlj(n2,lj2)
       nlj3 = SPB%nlj(n3,lj3)
       nlj4 = SPB%nlj(n4,lj4)
       DBD%TE2B(nlj1,nlj2,nlj3,nlj4,JJ) = DBD_TE
       DBD%TE2B(nlj3,nlj4,nlj1,nlj2,JJ) = DBD_TE
      goto 291
301   continue
      close(210)

      else
       stop ' Two-body DBD matrix elements are not initialized properly'
      endif
!   ..........................

      write(*,*) ' The TBMEs for NLDBD has been read !'
      return

      END


