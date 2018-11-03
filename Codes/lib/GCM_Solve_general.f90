      subroutine solution_general() 
      use VAPHFB_Par
      implicit none
      integer nmaxdj,maxmp,k1m,k2m
      integer jproj,k1,k2,kmax,iki0
      integer iqk1,iqk2,ikm,iq1,iq2
!    ......................................................... loop over
!    angular momentum J
      do jproj = GCM%Jmin,GCM%Jmax,1 !GCM%Jdf
         if(jproj.eq.1) cycle 

         write(*,102)
         write(*,201)
!        if(iv(jproj).gt.zero) then         ! for even spin
!           k1max  =  jproj
!        endif
!        if(iv(jproj).lt.zero) then       ! for odd spin
!           k1max  =  abs(jproj-1)
!        endif

!              do k1=-kmax,kmax,1   ! K-value
!              do k2=-kmax,kmax,1
!    ......................................................... print the
!    diagonal element 
                 call table_general(jproj,GCM%Iexst)
!              enddo
!              enddo
!    ......................................................... print the
!    norm overlap 
      print 400
      print 402
!     print 403
      write(*,102)
      maxmp = GCM%qkmax(jproj)
      do iqk1 = 1,maxmp   ! K-value
         iq1  = GCM%iq(iqk1,jproj)
         k1   = GCM%ik(iqk1,jproj)
         write(*,401) k1,iq1,(GCM%iexst(iqk1,iqk2),iqk2=1,iqk1)
      enddo
      !GCM%nmaxdi(iis)  = GCM%NOQ*(2*iis+1) 
      nmaxdj = GCM%nmaxdi(jproj)
!     ...........................................................
!     solution of HWG equation
      !print *,' presolve GCM: dim=',nmaxdj
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
  401 format(2x,2i3,1x,100i3)
  402 format ('    K  iq')
      
      return
      end
