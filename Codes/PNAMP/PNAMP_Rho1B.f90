!      
        subroutine ZRho1B_Angles(iangle,ZRO,Znorm,ZOBDJ0,NLEV)
        USE VAPHFB_PAR
        implicit none
        integer   mitf,iangle,NLEV 
        integer   k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 ZRO,ZOBDJ0,Znorm
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION ZOBDJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2

!     ............. print out the one-body density

        do k1=1,NLEV
        do k2=1,NLEV
           lj1 = tnljm%lj(k1)
           lj2 = tnljm%lj(k2)
           twoj1  = tnljm%twoj(k1)
           twoj2  = tnljm%twoj(k2)
           twom1  = tnljm%twom(k1)
           twom2  = tnljm%twom(k2)
           it1  = mitf(tnljm%t(k1))        ! -1 for p
           it2  = mitf(tnljm%t(k2))
           n1   = tnljm%n(k1)             ! n=1,2,3,...
           n2   = tnljm%n(k2)

           if(lj1 .ne. lj2 .or. twom1 .ne. twom2 .or. it1 .ne. it2) cycle
            ZOBDJ0(it1,lj1,n1-1,n2-1)= ZOBDJ0(it1,lj1,n1-1,n2-1) &
     &        + znorm*ZRO(k1,k2)/(dsqrt(twoj1+1.d0)*PNP%NFOM**2)  

!         if(abs(ZRO(k1,k2)) .gt. 1.d-4) &
!     &   write(*,'(4i4,2f12.5,5x,4f8.5)') it1,lj1,n1,n2,ZRO(k1,k2),znorm,ZOBDJ0(it1,lj1,n1-1,n2-1)
!    &   write(*,'(4i4,2f12.5,5x,2f12.5)') it1,lj1,n1,n2,ZRO(k1,k2),znorm
        enddo
        enddo
        end


        subroutine ZRho1B_Integral(iangle,ZOBDJ0,J,zff,ZRho1BJ0)
        USE VAPHFB_PAR
        implicit none
        integer   iangle
        integer   J,k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 zff,ZRho1BJ0,ZOBDJ0
        DIMENSION ZOBDJ0  (0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        DIMENSION ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2

       do it1=0,1
       do lj1=0,HO%ljmax
       do n1=0,HO%nMax
       do n2=0,HO%nMax
         ZRho1BJ0(it1,lj1,n1,n2)= ZRho1BJ0(it1,lj1,n1,n2) &
     &                        +weightJ(J)*zff*ZOBDJ0(it1,lj1,n1,n2) 

!     &                        +(2.d0*J+1.d0)*zff*ZOBDJ0(it1,lj1,n1,n2) &
!     &                           /(8.d0*pi**2)
!        if(abs(ZOBDJ0(it1,lj1,n1,n2)) .gt. 1.d-4) write(*,*) it1,lj1,n1,n2,ZOBDJ0(it1,lj1,n1,n2) 
!        if(abs(ZRho1BJ0(it1,lj1,n1,n2)) .gt. 1.d-4) &
!     &  write(*,*) it1,lj1,n1,n2,ZRho1BJ0(it1,lj1,n1,n2),zff 
       enddo
       enddo
       enddo
       enddo
       end

        subroutine Write_ME1B(ie,znorm,ZRho1BJ0)
        USE VAPHFB_PAR
        implicit none
        integer    ie,it,lj,n1,n2 
        integer    SPB_twoj_lj
        complex*16 znorm,ZRho1BJ0
        real*8     AE
        DIMENSION ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
!        real*8, parameter :: CHOP=1.d-20
!        real*8  eps
!        eps=3.d-14
        if (abs(znorm) .le. CHOP) then
           print *, 'warning: Znorm is too small'
           return 
        endif
 
        do it=0,1
        do lj=0,HO%ljmax
        do n1=0,HO%nMax
        do n2=0,HO%nMax
           AE= dreal(ZRho1BJ0(it,lj,n1,n2)/znorm)
           if(abs(AE) .gt. CHOP) write(ie,'(4i5,f20.10)') it,lj,n1,n2,AE/sqrt(SPB_twoj_lj(lj)+1.d0)
        enddo
        enddo
        enddo
        enddo
        write(*,'(a,i2)') '....The 1B-DME has been written to fort.',ie 
        end










       subroutine Write_Rho1B4IMSRG(OBD_Save)
       USE VAPHFB_Par
       implicit none
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,icount
       REAL*8, DIMENSION(0:1,1:HO%nMax,0:HO%LMax*2,0:1,1:HO%nMax,0:HO%LMax*2) :: OBD_Save
       real*8  OBD_t,Trace(0:1)
!   .............................................. 1B density matrix for the core state
     if(.NOT. ALLOCATED(SPB%CType)) ALLOCATE(SPB%CType(0:HO%ljmax,0:HO%nmax,0:1))
!  
      SPB%CType(0:HO%ljmax,0:HO%nmax,0:1) = -1
      OBD_Save(0:1,1:HO%nMax,0:HO%LMax*2,0:1,1:HO%nMax,0:HO%LMax*2) = zero 
!     .............................. loop over isopsin (it=0,1) 
      if(Num_Rho1B.eq.0) then
         icount=0
      do it=0, 1 ! n, p 
         icc = 0
      do e0=0, HO%eMax                 ! e0=0, ..., eMax
      do n1=1, e0/2+1               ! e0=2*(n1-1)+L
         L1=e0-2*(n1-1)
         L1=2*L1                    ! L1 is doubled  
      do J1=L1+1, abs(L1-1), -2     ! J1 is doubled
         LJ = (L1+J1-1)/2
         if(icc.lt.Nucl%ncore(it)) then
!     .............................. only for the core part
         OBD_Save(it,n1,LJ,it,n1,LJ)=1.d0  ! n1=1,2,3,
         icount = icount + 1         ! counting totally nonzero matrix elements 
         icc    = icc    + (J1+1)    ! (2j+1) counting total core nucleons
         SPB%CType(LJ,n1-1,it) = 1   ! n1-1 starts from 0
         endif
      enddo
      enddo
      enddo
      enddo
!  ......................................... counting the number of 1b density matrix elements
    open(59,file=File%Rho1B,status='old')
 10 read(59,*,end=20) it, LJ, n1, n2, OBD_t ! it=0 (n), 1 (p) 
    if(abs(OBD_t)>CHOP) icount=icount+1
    goto 10
 20 continue
    close(59)

    else
      icount = Num_Rho1B
    endif
!  ...................
      open(59,file=File%Rho1B,status='old')
      open(9,file='Rho1B.dat',status='unknown')
      write(9,31)
      write(9,41) icount
      write(9,11)


 30   read(59,91,end=40) it, LJ, n1, n2, OBD_t ! it=0 (n), 1 (p) 
        j1 = 2*int(LJ/2)+1  ! j1 is doubled 
        l1 =   LJ+1         ! l1 is doubled 
        j2 = j1
        l2 = l1
        OBD_Save(it,n1+1,LJ,it,n2+1,LJ)= OBD_t  ! j1 is doubled 
     goto 30
 40  continue
     close(59)

!     ............. print out the one-body density
      do it=0, 1 ! n, p 
         Trace(it) = zero
      do e0=0, HO%eMax                 ! e0=0, ..., eMax
      do n1=1, e0/2+1               ! e0=2*(n1-1)+L
      do n2=1, e0/2+1               ! e0=2*(n1-1)+L
         L1=e0-2*(n1-1)
         L1=2*L1                    ! L1 is doubled  
      do J1=L1+1, abs(L1-1), -2     ! J1 is doubled
         LJ = (L1+J1-1)/2
        if(n1.eq.n2) Trace(it) = Trace(it) + (J1+1.d0)*OBD_Save(it,n1,LJ,it,n2,LJ)
         if(abs(OBD_Save(it,n1,LJ,it,n2,LJ)).ge.CHOP) &
    &    write(9,91) it, LJ, n1-1, n2-1, OBD_Save(it,n1,LJ,it,n2,LJ)
      enddo
      enddo
      enddo
      enddo
         write(*,*) 'Trace[Rho1B]=',Trace(it)
      enddo
      close(9)
 11  format('       ')
 31  format('*** 1b DME normalized by 1/(2j+1) from GCM calculation ***')
 41  format(i5)
 91  format(4i5,2f15.10)
      end        


        subroutine Read_ME1B(ie,ZRho1BJ0)
        USE VAPHFB_PAR
        implicit none
        Complex*16  me
        integer     ie,it, LJ, n1, n2
        Complex*16  ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2

        open(ie,file='Rho1B.dat',status='old')
        read(ie,*)
        read(ie,*) Num_Rho1B
        read(ie,*)
11      read(ie,*,end=21) it, LJ, n1, n2, me

!        if(abs(me-one).lt.0.1) me = zone
!        if(abs(me).lt.0.1)     me = zzero
        ZRho1BJ0(it,LJ,n1,n2) = me   ! n=0,1,2

        goto 11

21      continue
!      ..............
        call Check_Trace_ZRho1B(ZRho1BJ0)

        end    

       subroutine Check_Trace_ZRho1B(ZRho1BJ0)
       USE VAPHFB_Par
       implicit none
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,nn,np
       Complex*16  ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
       Complex*16  Trace(0:1)
      do it=0, 1 ! n, p 
         Trace(it) = zero
      do e0=0, HO%eMax                 ! e0=0, ..., eMax
      do n1=1, e0/2+1               ! e0=2*(n1-1)+L
      do n2=1, e0/2+1               ! e0=2*(n1-1)+L
         L1=e0-2*(n1-1)
         L1=2*L1                    ! L1 is doubled  
      do J1=L1+1, abs(L1-1), -2     ! J1 is doubled
         LJ = (L1+J1-1)/2
        if(n1.eq.n2) Trace(it) = Trace(it) + (J1+1.d0)*ZRho1BJ0(it,LJ,n1-1,n2-1)  ! 
      enddo
      enddo
      enddo
      enddo
      enddo

        nn = Nucl%nucleon(0) 
        np = Nucl%nucleon(1) 
        if(abs(Trace(0)-nn).lt.1.d-3) then
           write(*,*) 'Trace[Rho1B](n) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho1B](n)=',Trace(0)
        endif
        if(abs(Trace(1)-np).lt.1.d-3) then
           write(*,*) 'Trace[Rho1B](p) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho1B](p)=',Trace(1)
        endif
      end    
