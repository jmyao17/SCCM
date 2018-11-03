      subroutine solution_triaxial() 
      use VAPHFB_Par
      implicit none
      integer nmaxdj,maxmp,k1m,k2m
      integer jproj,k1,k2,kmax,iki0
      integer iqk1,iqk2,ikm,iq1,iq2
!    ......................................................... loop over
!    angular momentum J
      do jproj = GCM%Jmin,GCM%Jmax,GCM%Jdf
         if(jproj.eq.1) cycle 
         if(AMP%i3DAMP.eq.0) then
              kmax = 0
         else 
               if(iv(jproj).gt.zero)  kmax = jproj   ! even
               if(iv(jproj).lt.zero)  kmax = abs(jproj-1)  ! odd
         endif
        ! write(*,102)
        ! write(*,'(a,i2)') ' Angular momentum J=',jproj
         write(*,102)
         write(*,201)
        
              do k1=-kmax,kmax,2   ! K-value
              do k2=-kmax,kmax,2
                 k1m = k1-iki0   ! mesh point
                 k2m = k2-iki0
!    ......................................................... print the
!    diagonal element 
                 call table(k1,k2,k1m,k2m,jproj,GCM%Iexst)
              enddo
              enddo
!    ......................................................... print the
!    norm overlap 
      print 400
      print 402
!     print 403
      write(*,102)
      maxmp = GCM%NOQ
      do k1=iki0,kmax,2   ! K-value
         k1m = k1-iki0    ! mesh point
      do iq1 = 1,maxmp
         iqk1=iq1+k1m/2*maxmp
         write(*,401) k1,iq1,(GCM%iexst(iqk1,iqk2),iqk2=1,iqk1)
      enddo
      enddo
      nmaxdj = GCM%nmaxdi(jproj)
!     ...........................................................
!     solution of HWG equation
      call PreSolve_HWG(jproj,nmaxdj,.false.)
      open(31,file=File%FF,status='unknown')
      print *,' start to solve GCM '
      call Solve_HWG(jproj,nmaxdj,.false.)
  600 enddo !jproj 

  102 format(/,1x,74('-'),/)
  201 format ('  J  K1 K2 iq     beta   gamma     qtot       P00     E &
      &      ',           '   n^J(q,q)    <N>     <Z>     <J^2>')

  400 format(/,' existing matrix elements',                         &
     &         ' -- value 99 means missing, otherwise ln(overlap)',/)
  401 format(2x,2i3,1x,35i3)
  402 format ('    K  iq')

      return
      end
