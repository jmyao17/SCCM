
      subroutine Read_Rho1B_J0(ie,Rho1BJ0,lpr)
!     .....................................
!     Density of ground state
!     .....................................
        USE VAPHFB_Par
        implicit none
        logical  lpr
        integer  ie,it, LJ, n1, n2
        real*8   cme,Rho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2

        Rho1BJ0=zero   
11      read(ie,*,end=21) it, LJ, n1, n2, cme
        Rho1BJ0(it,LJ,n1,n2) = cme   ! n=0,1,2
        goto 11
21      continue
!      ..............
       if(lpr) call Check_Trace_Rho1B(Rho1BJ0)
!      ..............
        end


       subroutine Add_Rho1B_J0(iq1,iq2,Rho1BJ0,lpr)
       USE VAPHFB_Par
       implicit none
       logical lpr
       integer n12,iq1,iq2
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,nn,np
       real*8  Rho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2

!      n12 = 2/(1 + iq1/iq2)  ! 1 or 2
      n12 = 1 + iq1/iq2  ! 1 or 2
      do it=0, 1 ! n, p 
      do e0=0, HO%eMax                 ! e0=0, ..., eMax
      do n1=1, e0/2+1               ! e0=2*(n1-1)+L
      do n2=1, e0/2+1               ! e0=2*(n1-1)+L
         L1=e0-2*(n1-1)
         L1=2*L1                    ! L1 is doubled  
      do J1=L1+1, abs(L1-1), -2     ! J1 is doubled
         LJ = (L1+J1-1)/2
         GCMDens%Rho1BSum(it,LJ,n1-1,n2-1) = &
      &                             GCMDens%Rho1BSum(it,LJ,n1-1,n2-1) &
      &                           + GCM%FqJk(iq1,0,1)*GCM%FqJk(iq2,0,1)       & 
      &                            *GCM%njkkqq(0,iq1,iq2)/n12         & 
      &                            *(Rho1BJ0(it,LJ,n1-1,n2-1)+Rho1BJ0(it,LJ,n2-1,n1-1)) 
      enddo
      enddo
      enddo
      enddo
      enddo

!     ............
       end


        subroutine Check_Rho1B(lpr)
        USE VAPHFB_PAR
        implicit none
        logical    lpr
        integer    ie,it,lj,n1,n2
        integer    SPB_twoj_lj
        real*8     cme

        write(*,*) ' ... Check the 1B Density Matrix  '
        open(59,file='For_GCM_Rho1BSum.dat',status='unknown')
        do it=0,1
        do lj=0,HO%ljmax
        do n1=0,HO%nMax
        do n2=0,HO%nMax
           cme= GCMDens%Rho1BSum(it,lj,n1,n2)
           if(abs(cme) .gt. CHOP) write(59,'(4i5,f12.8)') it,lj,n1,n2,cme
        enddo
        enddo
        enddo
        enddo
        close(59)
        if(lpr) call Check_Trace_Rho1B(GCMDens%Rho1BSum)
        end


       subroutine Check_Trace_Rho1B(Rho1BJ0)
       USE VAPHFB_Par
       implicit none
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,nn,np
       real*8  Rho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
       real*8  Trace(0:1)
      do it=0, 1 ! n, p 
         Trace(it) = zero
      do e0=0, HO%eMax                 ! e0=0, ..., eMax
      do n1=1, e0/2+1               ! e0=2*(n1-1)+L
      do n2=1, e0/2+1               ! e0=2*(n1-1)+L
         L1=e0-2*(n1-1)
         L1=2*L1                    ! L1 is doubled  
      do J1=L1+1, abs(L1-1), -2     ! J1 is doubled
         LJ = (L1+J1-1)/2
        if(n1.eq.n2) Trace(it) = Trace(it) + (J1+1.d0)*Rho1BJ0(it,LJ,n1-1,n2-1)  ! 
      enddo
      enddo
      enddo
      enddo
      enddo
        nn = Nucl%nucleon(0)
        np = Nucl%nucleon(1)
        if(abs(Trace(0)-nn).lt.1.d-3) then
           write(*,*) ' Trace[Rho1B](n) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho1B](n)=',Trace(0)
        endif
        if(abs(Trace(1)-np).lt.1.d-3) then
           write(*,*) ' Trace[Rho1B](p) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho1B](p)=',Trace(1)
        endif
      end


!     ........................................
       subroutine Write_Rho1B4IMSRG(OBD_Save)
!     ........................................
       USE VAPHFB_Par
       implicit none
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,icount
       REAL*8, DIMENSION(0:1,1:HO%nMax,0:HO%LMax*2,0:1,1:HO%nMax,0:HO%LMax*2) :: OBD_Save
       real*8  OBD_t,Trace(0:1)
!   .............................................. 1B density matrix for the core state

      if(.NOT. ALLOCATED(SPB%CType)) ALLOCATE(SPB%CType(0:HO%ljmax,0:HO%nmax,0:1))
!  
      SPB%CType= -1
      OBD_Save = zero
!     .............................. loop over isopsin (it=0,1) 
!     initial value of Num_Rho1B is 0  
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
      open(111,file='For_GCM_Rho1BSum.dat',status='unknown')
 10   read(111,*,end=20) it, LJ, n1, n2, OBD_t ! it=0 (n), 1 (p) 
      if(abs(OBD_t)>CHOP) icount=icount+1
      goto 10
 20   continue
      close(111)
!  ...................
      else
       icount = Num_Rho1B
      endif
!  ...................
      open(9,file='GCM_Rho1BSum.dat',status='unknown')
      write(9,31)
      write(9,41) icount
      write(9,11)
      open(112,file='For_GCM_Rho1BSum.dat',status='old')
 30   read(112,*,end=40) it, LJ, n1, n2, OBD_t ! it=0 (n), 1 (p) 
        j1 = 2*int(LJ/2)+1  ! j1 is doubled 
        l1 =   LJ+1         ! l1 is doubled 
        j2 = j1
        l2 = l1
        OBD_Save(it,n1+1,LJ,it,n2+1,LJ)= OBD_t  ! j1 is doubled 
      goto 30
 40  continue
     close(112)

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
         write(*,*) ' Trace[Rho1B]=',Trace(it)
      enddo
      close(9)
 11  format('       ')
 31  format('*** 1b DME normalized by 1/(2j+1) from GCM calculation ***')
 41  format(i5)
 91  format(4i5,2f15.10)
      end


        subroutine diag_Rho1B(lpr)
        USE VAPHFB_PAR
        implicit none
        logical    lpr
        integer    ie,it,lj,n1,n2,l2,j2
        integer    SPB_twoj_lj
        real*8     cme,x(1:HO%nMax+1,1:HO%nMax+1),aa(1:HO%nMax+1,1:HO%nMax+1)
        real*8     d(1:HO%nMax+1),e(1:HO%nMax+1)
        real*8     trace

        write(*,*) ' ... Diagonalizing the 1B Density Matrix  '
        open(69,file='GCM_Rho1BSum_natural.dat',status='unknown')
        do it=0,1
        trace = zero
        do lj=0,HO%ljmax
           l2 =   lj+1         ! l1 is doubled 
           j2 = 2*int(lj/2)+1  ! j1 is doubled 
        do n1=0,HO%nMax
        do n2=0,HO%nMax
           aa(n1+1,n2+1) = -GCMDens%Rho1BSum(it,lj,n1,n2)
        enddo
        enddo
          call sdiag(HO%nMax+1,HO%nMax+1,aa,d,x,e,1) 
          do n1=1,HO%nMax+1
             if(abs(d(n1)).gt.1.d-4) &
             write(69,'(2i6,a2,i3,a2,2f12.8)') &
                  it,n1-1,HO%tl(l2/2),j2,'/2',-d(n1),abs(d(n1))*(j2+1.0)
             trace = trace + abs(d(n1))*(j2+1.0)
          enddo
        enddo
         write(69,102) 
         write(69,'(a,f12.8)') 'Trace[Rho]=',trace
         write(69,102) 
        enddo
        close(69)
        write(*,*) ' The eigenvalues of 1B Density Matrix is stored into file: GCM_Rho1BSum_natural.dat'
 102    format(/,1x,74('-'),/)
        end 
