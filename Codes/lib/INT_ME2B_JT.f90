
!     read 2N in JT scheme from Heiko's file
! ..............................................................................
      SUBROUTINE INT_ME2J_JT()
! ..............................................................................
      USE VAPHFB_Par 
      IMPLICIT NONE
! ..............................................................................
!   This subroutine is to count the Two-Particle States
!   |(nlj)_1 (nlj)_2; J,PP,tt> in the block b=(J, PP,tt)
!   ...............
!   tt=0, 1, 2
!   PP=0, 1
!   e=2*n+L
!   e=1, 2, ...,eMax
!   L=0, 1, 2, Lmax
!   n=0, 1, 2, nmax=eMax/2
! ..............................................................................

      integer iphaseJT 
      real*8 p1p2,v2me,tot_me,d_ab,d_cd,AN_AB,AN_CD,delta_JT
      character(100) find_file,headline

      TYPE ME2J_i
        INTEGER iMax
        Real*8, dimension(:), allocatable :: me,tpp 
!        Real*8, dimension(:), allocatable :: nlj1,nlj2,nlj3,nlj4 
      END Type ME2J_i
      TYPE(ME2J_i)  pME2J
 
      INTEGER, PARAMETER :: nljMax0=100 
      TYPE ME2JH_nlj
        INTEGER e(0:nljMax0)
        INTEGER n(0:nljMax0)
        INTEGER lj(0:nljMax0)
        INTEGER l(0:nljMax0)
        INTEGER twoj(0:nljMax0)
        !INTEGER nlj(0:nMax,0:ljMax)
        !INTEGER idx(0:nljMax0,0:nljMax0,0:nljMax0,0:nljMax0)
      END Type ME2JH_nlj
      TYPE(ME2JH_nlj)  ME2JH
 
      INTEGER  nMax,ljMax,nljMax 
      INTEGER  nlj,e,n,L,lj,lMin  
      INTEGER  a, bb, c, bm, J, total,me
      INTEGER  t1, n1, l1, j1, lj1, id1 
      INTEGER  t2, n2, l2, j2, lj2, id2
      INTEGER  t3, n3, l3, j3, lj3, id3
      INTEGER  t4, n4, l4, j4, lj4, id4
      INTEGER  J12m, J12p
      INTEGER  LJ_1, LJ_2
      INTEGER  Maxdim, n_1234
      INTEGER  twoj,twojMin,twojMax,twol
      INTEGER  nlj1,nlj2,nlj3,nlj4
      INTEGER  J_min, J_max, TT, MT, JJ, PP
      INTEGER  i,nlj4max
      INTEGER  isum,iMax
      INTEGER  ia,ib,ic,id
      integer  endpoint

      character*100 ME_DIR,mefile

      integer iterm

      open(36,file=trim(INT_DIR)//File%f2bJTMT,status='unknown')
!     .............................determination of nljMax
      nljMax=-1 
      nlj   =0 !-1 
      !ME2JH%idx = -1
      do e=0,HO%eMax
         lMin = mod(e,2)            ! 0 (even eMax); or 1 (odd eMax)
         do l=lMin,min(e,HO%lMax),2  ! L is not doubled
            n=(e-l)/2   ! 0,1,..
            if(n.gt.HO%nMax) cycle
            twol = 2*l
            twojMin = abs(twol-1) 
            twojMax = twol+1 
            do twoj=twojMin, twojMax, 2
               lj = (twol + twoj -1 )/2
               nlj=nlj + 1 
               if(nlj.gt.nljMax0) &
     &         stop '*Error*: nljMax0 should be larger !!'
                ME2JH%e(nlj)     = e
                ME2JH%n(nlj)     = n
                ME2JH%lj(nlj)    = lj
                ME2JH%l(nlj)     = l
                ME2JH%twoj(nlj)  = twoj
                !ME2JH%nlj(n,lj) = nlj
            enddo
          enddo
        enddo
      nljMax=nlj 
      write(*,*) ' nljMax=',nljMax
!     .............................. build up the index 
       i = -1 
       do nlj1=1,nljMax
       do nlj2=1,nlj1 ! nljMax
          if(ME2JH%n(nlj1)+ME2JH%n(nlj2).gt.HO%eMax*2) cycle

!      .........................
       do nlj3=1,nlj1 ! nljMax
!         loop  nlj4 = 0...nlj2 or nlj3 nlj3 >= nlj4 and nlj1*nljDim+nlj2 >= nlj3*nljDim+nlj4 required
!        (nlj3==nlj1 ? nlj2 : nlj3)
          nlj4max=nlj3
          if(nlj3.eq.nlj1) nlj4max=nlj2

       do nlj4=1,nlj4max

          if(ME2JH%n(nlj3)+ME2JH%n(nlj4).gt.HO%eMax*2) cycle

          !ME2JH%idx(nlj1,nlj2,nlj3,nlj4) = -1

          if(iv(ME2JH%l(nlj1)+ME2JH%l(nlj2))          &
     &    .ne.iv(ME2JH%l(nlj3)+ME2JH%l(nlj4))) cycle
!           if(mod(ME2JH%l(nlj1)+ME2JH%l(nlj2)          &
!     &       +ME2JH%l(nlj3)+ME2JH%l(nlj4),2).ne.0) cycle

           J_min=max(abs(ME2JH%twoj(nlj1)-ME2JH%twoj(nlj2))/2, &
     &              abs(ME2JH%twoj(nlj3)-ME2JH%twoj(nlj4))/2) 
           J_max=min(abs(ME2JH%twoj(nlj1)+ME2JH%twoj(nlj2))/2, &
     &              abs(ME2JH%twoj(nlj3)+ME2JH%twoj(nlj4))/2) 

          if(J_min.gt.J_max) cycle
          !ME2JH%idx(nlj1,nlj2,nlj3,nlj4) = i+1
!    ................. 
           do JJ=J_min,J_max          
           do TT=0,1,+1          
           do MT=-TT,TT          
              i = i + 1

          enddo
          enddo
          enddo
!    ................. 
          enddo
          enddo
          enddo
          enddo
       pME2J%iMax=i
       
       write(*,*) ' pME2J%iMax=',i
 20    format(4i3,i8)


!     .............. Read ME2J

      allocate(pME2J%me(0:pME2J%iMax))
      allocate(pME2J%tpp(0:pME2J%iMax))

      pME2J%me(:) = 0.d0 
      pME2J%tpp(:) = 0.d0 
!   ............ two parts: center of mass correction (2b) + original two-body

      do iterm =1 ,2

      if(iterm.eq.1) mefile=File%NN !'chi2b_srg0953_eMax06_hwHO020.me2j'
      if(iterm.eq.2) mefile=File%p1p2 !'tpp_eMax06.me2j'
      ME_DIR = find_file("ME_FILES",mefile)
      open(9,file=trim(ME_DIR)//trim(mefile),status='old')
      ! print *,ME_DIR
      read(9,*) headline
      !print *, headline
      isum = 0
      do i = 0,pME2J%iMax,10
         if(i+10 .le.pME2J%iMax) endpoint = 10 
         if(i+10 .gt.pME2J%iMax) endpoint = pME2J%iMax-i+1 
         if(iterm.eq.1) read(9,'(f12.7,9f13.7)') (pME2J%me(isum+j), j=1,endpoint)
         if(iterm.eq.2) read(9,'(f12.7,9f13.7)') (pME2J%tpp(isum+j), j=1,endpoint)
         isum = isum +10
       enddo

     ! ...........
      enddo ! iterms

       i = 0
       do nlj1=1,nljMax
            lj1 = ME2JH%lj(nlj1)   
            n1  = ME2JH%n(nlj1)
            ia = SPB%nlj(n1,lj1)

       do nlj2=1,nlj1 ! nljMax

            lj2 = ME2JH%lj(nlj2)   
            n2  = ME2JH%n(nlj2)
            ib = SPB%nlj(n2,lj2)

                     d_ab = zero
        if(nlj1.eq.nlj2) d_ab = one
!      .........................
          if(ME2JH%n(nlj1)+ME2JH%n(nlj2).gt.HO%eMax*2) cycle
          do nlj3=1,nlj1 ! nljMax
             nlj4max=nlj3
             if(nlj3.eq.nlj1) nlj4max=nlj2

            lj3 = ME2JH%lj(nlj3)   
            n3  = ME2JH%n(nlj3)
             ic = SPB%nlj(n3,lj3)

          do nlj4=1,nlj4max
 
             if(ME2JH%n(nlj3)+ME2JH%n(nlj4).gt.HO%eMax*2) cycle
            lj4 = ME2JH%lj(nlj4)   
            n4  = ME2JH%n(nlj4)
             id = SPB%nlj(n4,lj4)

                     d_cd = zero
            if(nlj3.eq.nlj4) d_cd = one

!             ME2JH%idx(nlj1,nlj2,nlj3,nlj4) = -1

          if(iv(ME2JH%l(nlj1)+ME2JH%l(nlj2))          &
     &    .ne.iv(ME2JH%l(nlj3)+ME2JH%l(nlj4))) cycle

           J_min=max(abs(ME2JH%twoj(nlj1)-ME2JH%twoj(nlj2))/2, &
     &              abs(ME2JH%twoj(nlj3)-ME2JH%twoj(nlj4))/2)
           J_max=min(abs(ME2JH%twoj(nlj1)+ME2JH%twoj(nlj2))/2, &
     &              abs(ME2JH%twoj(nlj3)+ME2JH%twoj(nlj4))/2)

          if(J_min.gt.J_max) cycle
!    ................. 
           do JJ=J_min,J_max
           do TT=0,1,+1
           do MT=-TT,TT   


              i = i + 1

            iphaseJT  = iv(JJ+TT)
            AN_AB     = sqrt((1.d0-d_ab*iphaseJT))/(1.d0+d_ab)
            AN_CD     = sqrt((1.d0-d_cd*iphaseJT))/(1.d0+d_cd)
            delta_JT  = AN_AB*AN_CD

         !   lj1 = ME2JH%lj(nlj1)    
         !   n1  = ME2JH%n(nlj1)    

! attention: in ME2J file, MT=-1 corresponding to nn. we convert it to +1 for nn here 

            p1p2   = delta_JT*2.d0*pME2J%tpp(i)/(Nucl%nucleon(2)*HO%b_osc**2)
            v2me   = delta_JT*pME2J%me(i)
            tot_me = v2me + p1p2

            if(abs(tot_me) .gt. CHOP)  &
            write(36,'(7i5,f15.8)') ia,ib,ic,id,JJ,TT,-MT,tot_me

          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo

          print *, ' reading ME2J is completed ......'

        close(36)
      return

      END 


