      subroutine Generate_ME1B4Q2(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
!     ...........................................
!     read Q2mu(1b) matrix elements (un-normal ordered)
!          from the files by Nathan. 
!     ...........................................
        USE VAPHFB_PAR
        implicit none
        real*8 Q2_2t,Q2_1t,Q20t,Q21t,Q22t
        DIMENSION Q2_2t(HO%NLEV,HO%NLEV)
        DIMENSION Q2_1t(HO%NLEV,HO%NLEV)
        DIMENSION Q20t(HO%NLEV,HO%NLEV)
        DIMENSION Q21t(HO%NLEV,HO%NLEV)
        DIMENSION Q22t(HO%NLEV,HO%NLEV)
        Integer t1,LJ,it,nlj,tz
        Integer n1,l1,j1,tz1
        integer n2,l2,j2,tz2,LJ1,LJ2
        real*8 fac,cg1,cg2,cg3,cg4,cg5,f1bJ(-1:1,1:HO%nmax,1:HO%nmax,0:HO%ljmax,0:HO%ljmax)
        Integer ii,jj,lk
        Real*8 ee_t,f1bm(HO%NLEV,HO%NLEV),f0b
        integer kk,ll,tk,nk,twojk,twomk,ljk,tl,nl,twojl,twoml,ljl
        character(len=100) :: file_nathan
        REAL(DP) :: CG
!     ............. initialization
      f1bJ(-1:1,1:HO%nmax,1:HO%nmax,0:HO%ljmax,0:HO%ljmax) = zero

!       file_nathan="../Int/Nathan_IMSRG_"//Input%cFlow//"_Q2init_"//Input%cIntID//&
      file_nathan="../Int/Nathan_IMSRG_"//Input%cFlow//"_Q2_"//Input%cIntID//&
          & "_"//Input%cValID//'_'//Input%chwHO//".me1b.dat"

       open(192,file='Nathan_Q20_'//Input%chwHO//'_'//&
      &    Input%cValID//'_me1b.dat',status='unknown')

       write(*,*) '    Read S.P.E. from',file_nathan

       open(21,file=file_nathan,status="old")
!     .................. read sp part of Q2 
       read(21,*)     
       read(21,*) f0b 
       read(21,*)     
       read(21,*)     
       read(21,*)    
!    ....................................................... 
!    n1  l1  2*j1  2*tz1     n2  l2  2*j2  2*tz2      T(1,2)
!    ....................................................... 
!    0    0    1   -1       0    2    3   -1           2.7618696
!    0    0    1   -1       0    2    5   -1           3.3567432
!    0    0    1   -1       1    2    3   -1          -0.4480555
!    0    0    1   -1       1    2    5   -1          -0.5485826
!    0    0    1   -1       2    2    3   -1           0.0210444
!    ....
!    tz=+1(n); -1(p)
!    ....................................................... 
 70    read(21,*,end=80) n1,l1,j1,tz1,n2,l2,j2,tz2, ee_t  ! n1 =0,1,2
          if(tz1.ne.tz2) stop ' Error: n-p mixing in 1b matrix elements'
          LJ1 = (2*l1 + j1 -1)/2
          LJ2 = (2*l2 + j2 -1)/2
          f1bJ(tz1,n1+1,n2+1,LJ1,LJ2) = ee_t
          f1bJ(tz1,n2+1,n1+1,LJ2,LJ1) = ee_t*iv(abs(j1-j2)/2)
!          print *, n1,l1,j1,tz1,n2,l2,j2,tz2, ee_t
          goto 70
 80    continue


!      .initialization
        DO kk=1,HO%NLEV
        DO ll=1,HO%NLEV
           Q2_2t(ll,kk)= 0.d0 
           Q2_1t(ll,kk)= 0.d0 
           Q20t(ll,kk)=  0.d0 
           Q21t(ll,kk)=  0.d0 
           Q22t(ll,kk)=  0.d0 
        enddo
        enddo
       
        DO ll=1,HO%NLEV
          tl = tnljm%t(ll)  ! -1(p); +1 (n) 
          nl = tnljm%n(ll)  ! 1,2,3,...
          twojl = tnljm%twoj(ll)
          twoml = tnljm%twom(ll)
          ljl   = tnljm%lj(ll)
          if(tl.ne.-1) cycle  ! only protons have nonzero matrix elements
!          write(500,'(6i5)') ll,tl,nl,ljl,twojl,twoml
        DO kk=1,HO%NLEV
          tk = tnljm%t(kk)  ! -1(p); +1 (n) 
          nk = tnljm%n(kk)  ! 1,2,3,...
          twojk = tnljm%twoj(kk)
          twomk = tnljm%twom(kk)
          ljk   = tnljm%lj(kk)

          if(tk.ne.tl) cycle
         
          fac = iv(abs(twojk-twomk)/2)/sq(5)

           cg1=fac*CG(twojl,twoml,twojk,-twomk,2*2,-2*2) 
           cg2=fac*CG(twojl,twoml,twojk,-twomk,2*2,-1*2) 
           cg3=fac*CG(twojl,twoml,twojk,-twomk,2*2,0) 
           cg4=fac*CG(twojl,twoml,twojk,-twomk,2*2,1*2) 
           cg5=fac*CG(twojl,twoml,twojk,-twomk,2*2,2*2)            

           Q2_2t(ll,kk)= f1bJ(tk,nl,nk,ljl,ljk)*cg1
           Q2_1t(ll,kk)= f1bJ(tk,nl,nk,ljl,ljk)*cg2
           Q20t(ll,kk)= f1bJ(tk,nl,nk,ljl,ljk)*cg3
           Q21t(ll,kk)= f1bJ(tk,nl,nk,ljl,ljk)*cg4
           Q22t(ll,kk)= f1bJ(tk,nl,nk,ljl,ljk)*cg5
        END DO
        END DO
        print *, 'Read me1b from Nathans file'

        return
        DO ll=1,HO%NLEV
        DO kk=1,HO%NLEV
          if(abs(Q20t(ll,kk)).gt.1.d-5) &
       &  write(192,'(2i4,f10.5)') ll,kk,Q20t(ll,kk) !f1bJ(tk,nl,nk,ljl,ljk)
        enddo
        enddo
        end subroutine

