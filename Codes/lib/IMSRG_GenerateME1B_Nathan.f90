      subroutine Generate_ME1B(f0b,f1bm)
!     ...........................................
!     read 1b matrix elements (un-normal ordered)
!          from the files by Nathan. 
!     ...........................................
        USE VAPHFB_PAR
        implicit none
        Integer t1,LJ,it,nlj,tz
        Integer n1,l1,j1,tz1
        integer n2,l2,j2,tz2
        real*8 Ascale,f1bJ(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax)
        Integer ii,jj,lk
        Real*8 ee_t,f1bm(HO%NLEV,HO%NLEV),f0b
        integer kk,ll,tk,nk,twojk,twomk,ljk,tl,nl,twojl,twoml,ljl
        character*14 ctemp

!     ............. initialization
      f1bJ(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax) = zero

      open(41,file='chi2b3b.int',status='old')
      read(41,'(a14,a)') ctemp,File%IMSRG_Hme1b  
      read(41,'(a14,a)') ctemp,File%IMSRG_Hme2b
      close(41)  
      print *,File%IMSRG_Hme1b,File%IMSRG_Hme2b

!      file_nathan="../Int/Nathan_IMSRG_"//Input%cFlow//"_SPE_"//Input%cIntID//&
!         & "_"//Input%cValID//'_'//Input%chwHO//".dat"
!
       write(*,*) '    Read S.P.E. from',File%IMSRG_Hme1b
       open(21,file=File%IMSRG_Hme1b,status="old")
!     .................. read sp part of H
       read(21,*)     
       read(21,*) f0b !H%E0
       read(21,*)     
       read(21,*)     
       read(21,*)    
!    ....................................................... 
!    n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      T(1,2)
!    ....................................................... 
!    tz=+1(n); -1(p)
!    ....................................................... 
 70    read(21,*,end=80) n1,l1,j1,tz1,n2,l2,j2,tz2, ee_t
          if(tz1.ne.tz2) stop ' Error: n-p mixing in 1b matrix elements'
          LJ = (2*l1 + j1 -1)/2
          if(LJ.ne.(2*l2 + j2 -1)/2) &
        & stop ' Error: LJ number in 1b matrix elements' 
          
          f1bJ(tz1,n1+1,n2+1,LJ) = ee_t
          f1bJ(tz1,n2+1,n1+1,LJ) = ee_t
!          print *, n1,l1,j1,tz1,n2,l2,j2,tz2, ee_t
          goto 70
 80    continue

!      .......... print out the ME1B
!  //*** loop over quantum numbers
!  for (t=0; t<=1; t++) {
!        for (lj=0; lj<=SPB.ljMax; lj++) {
!          for (n=0; n<=SPB.nMax_lj[lj]; n++) {
!                for (nn=0; nn<=SPB.nMax_lj[lj]; nn++) {
!                  me = ME1B->me[t][lj][n][nn];
       
       print *, ' print the ME1B into',File%f1b
       open(22,file=File%f1b,status="unknown")
       write(22,*) H%E0
       do it=0,1
          tz=iv(it)
       do nlj = 1, HO%nljmax
          n1   = SPB%n(nlj)  ! 0,1,2,..
          LJ  = SPB%lj(nlj)
          L1  = Int((LJ+1)/2)
       do n2 = 0, HO%nmax-1
          if(2*n2+L1.gt.HO%emax) cycle
          if(abs(f1bJ(tz,n1+1,n2+1,LJ)).gt.1.d-10) &
        & write(22,'(4i8,f20.8)') it, LJ, n1, n2, f1bJ(tz,n1+1,n2+1,LJ)
         if(abs(f1bJ(tz,n1+1,n2+1,LJ)).gt.1000) write(*,'(4i8,f15.8)') it, LJ, n1, n2, f1bJ(tz,n1+1,n2+1,LJ)
       enddo
       enddo
       enddo

        !initialization
        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            f1bm(jj,ii)=zero
         end do
        end do
        DO kk=1,HO%NLEV
        DO ll=1,HO%NLEV

          tk = tnljm%t(kk)  ! -1(p); +1 (n) 
          nk = tnljm%n(kk)  ! 1,2,3,...
          twojk = tnljm%twoj(kk)
          twomk = tnljm%twom(kk)
          ljk   = tnljm%lj(kk)

          tl = tnljm%t(ll)  ! -1(p); +1 (n) 
          nl = tnljm%n(ll)  ! 1,2,3,...
          twojl = tnljm%twoj(ll)
          twoml = tnljm%twom(ll)
          ljl   = tnljm%lj(ll)
          if(tk .ne. tl .or. ljk.ne.ljl .or. twomk .ne. twoml) cycle
          f1bm(ll,kk) = f1bJ(tk,nl,nk,ljk)
!          write(*,*) tk,ll,kk,f1bm(ll,kk)
        END DO
        END DO
        end subroutine

