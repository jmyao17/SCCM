!______________________________________________________________________________ 
      subroutine PreSolve_HWG(jtot,nmaxd,lpr)

!.............................................................................. 
!           Solution of Hill-Wheel equation 
!       The standard method is adopted to solve
!       the HW equation. More details can be found
!       in <The nuclear many-body problem> by P.Ring & P.Schuck
!.............................................................................. 
      use VAPHFB_Par
      implicit real*8 (a-h,o-z)
      real*8 NN,HH
      logical lpr
!---------------------------------------------------------------------
      real*8    checkzero
      dimension HH(nmaxd,nmaxd),NN(nmaxd,nmaxd),GG(nmaxd,nmaxd),  &
     &          DD(nmaxd,nmaxd),FF(nmaxd,nmaxd)
      DIMENSION WW(nmaxd,nmaxd),RR(nmaxd,nmaxd)
      dimension E(nmaxd),EN(nmaxd),Z(nmaxd)

      dimension CheckE(nmaxd,nmaxd)
   10 format(3x,5f18.12)
  100 format ('cden(',2i4,')=',2f12.8)
  101 format(3i4,f12.5,5f9.5)
  102 format(/,1x,74('-'),/)
  200 format (//,' eigenvalues of the norm ',/)
  201 format (1h ,1p7e11.3)
  202 format (/,' cutoff value of the norm eigenvalues ',1pe10.2)
  203 format (//,' eigenvalues of the hamiltonian ',/)
  401 format(2x,2i3,1x,3f8.3,10f10.5)
  402 format ('    K  iq   beta  gamma    P00  '   &
     &   '      g1       g2       g3 ')
  403 format (//,'  eigenvectors ',/)
  404 format ('    K  iq   beta  gamma    P00  '   &
     &   '     fn1      fn2      fn3 ')

!     .............................................. initialization  
      if(GCM%kmax .gt.10) GCM%kmax = 10
      N = nmaxd
      NN(:,:)   = zero
      HH(:,:)   = zero
      checkzero = zero
      do iq1=1,N
         do iq2=1,N   
            !if(abs(GCM%njkkqq(jtot,iq2,iq1)) .lt. 1.d-4) cycle 
            NN(iq2,iq1) = GCM%njkkqq(jtot,iq2,iq1)
            HH(iq2,iq1) = GCM%hjkkqq(jtot,iq2,iq1) !*NN(iq2,iq1) !GCM%njkkqq(jtot,iq2,iq1)
            !write(*,'(2i5,2f12.8)') iq2,iq1,NN(iq2,iq1),HH(iq2,iq1)

            if(abs(NN(iq2,iq1)).gt.checkzero) &
            checkzero = NN(iq2,iq1)
         enddo !q2
      enddo !q1 
 
      if(abs(checkzero).lt.1.d-6) then
        print *, ' THIS STATE DOES NOT EXIST !'
        return ! this state does not exist ...
      endif 
!-----------------------------------------------------------------
      GCM%variance(:,:)=100
      do NOS = 1, N
         !call HNDIAG(NMAXD,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,Z,NOS,0,IFL)
         call HNDIAG(NMAXD,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,Z,NOS,3,IFL)
       do k=1,min(GCM%kmax,NOS) 
           CheckE(k,NOS) = E(k)
          do iqk=1,N
             GCM%FqJk(iqk,jtot,k) = FF(iqk,k)
          enddo ! iqk 
             GCM%EJk(jtot,k)     = E(k)
             call mean_variance(jtot,k,q1,q2)
             if(abs(GCM%variance(jtot,k)).gt.abs(q2)) then
                  GCM%NOS_J(k,jtot) = NOS
                  GCM%variance(jtot,k) = q2
                  GCM%mean(jtot,k)     = q1
             endif
             !if(jtot.eq.4) print *, jtot,k,GCM%NOS_J(k,jtot)
       enddo ! k
      enddo ! NOS
!   ................... print out the energy for checking convergence
      write(*,102)
      write(*,*) ' Convergence Check: E vs No. of Natural States '
      write(*,*) ' NOS    E_1       E_2     ...'
      jout=jtot+200
      do NOS=1,N
         write(*,'(i5,100f10.3)') NOS, (CheckE(k,NOS), k=1,NOS)
         write(jout,'(i5,100f10.3)') NOS, (CheckE(k,NOS), k=1,NOS)
      enddo ! NOS 

      write(*,102)
      if(Input%NOSJ(jtot).ge.1)  GCM%NOS_J(1,jtot)=Input%NOSJ(jtot)
      GCM%NOS_J(2,0) = 8
      do k=1, min(GCM%kmax,N)
         write(*,'(a,2i3,a,i3)') &
     &   '  Number of configurations for the state (J,k)=',jtot,k,&
     &   ' is ',GCM%NOS_J(k,jtot)
!     ..................................
      enddo ! k 
   
      return
      end
!---------------------------------------------------------------------
      subroutine Solve_HWG(jtot,nmaxd,lpr)
!......................................... 
!           Solution of Hill-Wheel equation 
!       The standard method is adopted to solve
!       the HW equation. More details can be found
!       in <The nuclear many-body problem> by P.Ring & P.Schuck
!.............................................................................. 
      use VAPHFB_Par
      implicit real*8 (a-h,o-z)
      real*8 NN,HH
      logical lpr
      dimension HH(nmaxd,nmaxd),NN(nmaxd,nmaxd),GG(nmaxd,nmaxd),  &
     &          DD(nmaxd,nmaxd),FF(nmaxd,nmaxd)
      DIMENSION WW(nmaxd,nmaxd),RR(nmaxd,nmaxd)
      dimension E(nmaxd),EN(nmaxd),Z(nmaxd)
   10 format(3x,5f18.12)
  100 format ('cden(',2i4,')=',2f12.8)
  101 format(3i4,f12.5,5f9.5)
  102 format(/,1x,74('-'),/)
  200 format (//,' eigenvalues of the norm ',/)
  201 format (1h ,1p7e11.3)
  202 format (/,' cutoff value of the norm eigenvalues ',1pe10.2)
  203 format (//,' eigenvalues of the hamiltonian ',/)
  401 format(2x,2i3,1x,3f8.3,10f10.5)
  402 format ('    K  iq   beta  gamma    P00  '   &
     &   '     g1      g2      g3 ')
  404 format ('    K  iq   beta  gamma    P00  '   &
     &   '     f1      f2      f3 ')

  403 format (//,'  eigenvectors ',/)
!    ............... re-initialization of NN and HH
      N=nmaxd
      do iq1=1,N
         do iq2=1,N
            NN(iq2,iq1) = GCM%njkkqq(jtot,iq2,iq1)
            HH(iq2,iq1) = GCM%hjkkqq(jtot,iq2,iq1) !*GCM%njkkqq(jtot,iq2,iq1)
         enddo !q2
!            if(jtot.eq.0) write(*,'(100f10.3)') (HH(iq2,iq1),iq2=1,N) 
      enddo !q1 
 
!     .............. diag. H
      write(*,102)
      do k=1,min(GCM%kmax,N)
         NOS = GCM%NOS_J(k,jtot)
         call HNDIAG(NMAXD,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,Z,NOS,3,IFL)

!     ...........................
!     ATTENTION: avoid over-write by introducing FqJk and GJkq
!     ...........................
          do iqk=1,N
             GCM%GJkq(iqk,jtot,k) = RR(iqk,k)
             GCM%FqJk(iqk,jtot,k) = FF(iqk,k)
             !if(jtot.eq.0) write(*,'(3i4,f12.8)') iqk,jtot,k,GCM%FqJk(iqk,jtot,k) 
          enddo ! iqk 
             GCM%EJk(jtot,k)     = E(k)
             write(*,'(a,i4,a,2f10.3)') &
             ' With NOS=',NOS,' Energy of States:', GCM%EJk(jtot,k)
      enddo ! k

!     ............................................  RR: g(q)
      write(*,102)
      print 402
      do iqk1=1,N

         if(Input%icr.eq.0) then
         ! no cranking
            k1 =2*Int((iqk1-1)/GCM%NOQ)          ! k value
            if(iv(jtot).lt.zero) k1 = k1 + 2   ! if jtot = 3, 5, 7, ...; k value
            iq1 = mod(iqk1,GCM%NOQ)              ! mesh point only in q-space
            if(iq1.eq.0) iq1 = GCM%NOQ 
          else
         ! with cranking
            k1  = GCM%ik(iqk1,jtot)
            iq1 = GCM%iq(iqk1,jtot)
          endif

         write(*,401) k1,iq1,Const%beta2t_mesh(iq1),Const%gamma2t_mesh(iq1), &
     &   Const%P00_mesh(iq1),(GCM%GJkq(iqk1,jtot,ik),ik=1,min(N,GCM%kmax))
      enddo ! iqk1     
!-----------------------------------------------------------------
      print 404
        do iqk1=1,N

         if(Input%icr.eq.0) then
            k1 =2*Int((iqk1-1)/GCM%NOQ)          ! k value
            if(iv(jtot).lt.zero) k1 = k1 + 2   ! if jtot = 3, 5, 7, ...;k value
            iq1 = mod(iqk1,GCM%NOQ)              ! mesh point only in q-space
            if(iq1.eq.0) iq1 = GCM%NOQ
          else
            k1  = GCM%ik(iqk1,jtot)
            iq1 = GCM%iq(iqk1,jtot)
          endif

         write(*,401) k1,iq1,Const%beta2t_mesh(iq1),Const%gamma2t_mesh(iq1), &
     &              Const%P00_mesh(iq1),(GCM%FqJk(iqk1,jtot,ik),ik=1,min(N,GCM%kmax))
!      ....... print out the F for the ground state
          if(jtot.eq.0) write(31,'(i3,2f15.8)') iqk1,GCM%FqJk(iqk1,jtot,input%J0k),NN(iqk1,iqk1)
         enddo ! iqk1 
          if(jtot.eq.0) print*, ' The FF for the G.S. is written into',File%FF
      return
      end


      subroutine mean_variance(ji,li,q1,q2)
!    ......
!     ji for spin
!     li =1, 2, ...
!    .....
      USE VAPHFB_Par
      implicit real*8(a-h,o-z)

      real*8  aa(GCM%NOQ*(JJmax*2+1),GCM%NOQ*(JJmax*2+1))
      real*8  bb(GCM%NOQ*(JJmax*2+1),GCM%NOQ*(JJmax*2+1))
      real*8  s(GCM%NOQ*(JJmax*2+1),GCM%kmax,0:JJmax)
      real*8  y(GCM%NOQ*(JJmax*2+1)),zz(GCM%NOQ*(JJmax*2+1)),w,x

        im   = GCM%nmaxdi(ji)
        ek = GCM%EJk(ji,li)
        do iqk1 = 1, im
           s(iqk1,li,ji)  = GCM%FqJk(iqk1,ji,li)
        enddo !  iqk1

        do iqk1 = 1, im
        do iqk2 = 1, im
           aa(iqk1,iqk2)  = GCM%hjkkqq(ji,iqk1,iqk2) !*GCM%njkkqq(ji,iqk1,iqk2)
           bb(iqk1,iqk2)  = GCM%njkkqq(ji,iqk1,iqk2)
           !print *, iqk1,iqk2,aa(iqk1,iqk2)-ek*bb(iqk1,iqk2)
        enddo !  iqk1
        enddo !  iqk2
        do i=1,im
           y (i) = 0.0d0
           zz(i) = 0.0d0
           do j=1,im
              x = 0.0
!     ...........................
            do l=1,im
               x = x + (aa(j,l)-ek*bb(j,l)) * s(l,li,ji) ! x=  sum_l (H_jl - E N_jl) f_l
            enddo  ! l
!     ...........................
               w    = aa(i,j)-ek*bb(i,j)
               y(i) = y(i) + w*x               ! y_i = sum_j (H_ij - E N_ij) * x
               zz(i)= zz(i)+ w*s(j,li,ji)
          enddo ! j
        enddo  ! i
!     .............................
        q1 = 0.0d0
        q2 = 0.0d0
        do i=1,im
          q1 = q1 + y (i)*y(i)      
          q2 = q2 + zz(i)*zz(i)       ! sum_i [ sum_j (H_ij - E N_ij)* f_j ]**2
        enddo ! i 
        end

