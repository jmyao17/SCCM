!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine with the single particle energies expressed in the uncoupled basis n,l,m_l,m_s_mt!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine kinetic(Esp,Ascale)
        USE VAPHFB_PAR
        use omp_lib
        implicit none 
        Integer t1,LJ
        real*8 Ascale,ee(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax)
        Integer ii,jj,lk 
        Real*8 ee_t,Esp(HO%NLEV,HO%NLEV)
        integer kk,ll,tk,nk,twomk,ljk,tl,nl,twoml,ljl
        integer nlj1,nlj2,lj1,lj2,n1,n2,l1,l2
         
       if(Input%IntJT.eq.3) then
!     ............. computing the s.p.e  
           do nlj1 = 1, HO%nljmax
              n1 = SPB%n(nlj1)
              l1 = SPB%l(nlj1)
              lj1= SPB%lj(nlj1)

           do nlj2 = 1, HO%nljmax
              n2 = SPB%n(nlj2)
              l2 = SPB%l(nlj2)
              lj2= SPB%lj(nlj2)

              if(lj1 .ne. lj2) cycle
              ee_t =0.d0
              if(n1.eq.n2)   ee_t = HO%hb0*(2*n1+l1+1.5)/HO%b_osc**2 
              if(n1.eq.n2-1) ee_t = HO%hb0*sqrt((n1+1)*(n1+l1+1.5))/HO%b_osc**2 
              if(n1.eq.n2+1) ee_t = HO%hb0*sqrt((n2+1)*(n2+l2+1.5))/HO%b_osc**2 

              ee(1,n1+1,n2+1,lj1)  = ee_t*(one-one/Nucl%nucleon(2))
              ee(-1,n1+1,n2+1,lj1) = ee_t*(one-one/Nucl%nucleon(2))

               ! write(*,*) lj1, n1, n2, ee(1,n1+1,n2+1,lj1)
             enddo
          enddo


       else
!     ............. read the s.p.e from file 
           open(21,file=trim(INT_DIR)//File%IMSRG_Hme1b,status="old") 
!     ............. initialization
           ee(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax) = zero
           write(*,*) '    Read S.P.E. from',trim(INT_DIR)//File%IMSRG_Hme1b
!     .................. read sp part of H
!          t1=0 (n); 1 (p); 2 (n+p)
           read(21,*) H%E0
 70        read(21,*,end=80) t1, LJ, n1, n2, ee_t
           if(Input%IntType .eq. 0) ee(iv(t1),n1+1,n2+1,LJ) = ee_t              
           if(Input%IntType .eq. 1 .and. Input%iCOM.eq.1) ee(iv(t1),n1+1,n2+1,LJ) = ee_t*(one-one/Nucl%nucleon(2))
           if(Input%IntType .eq. 1 .and. Input%iCOM.eq.2) ee(iv(t1),n1+1,n2+1,LJ) = ee_t
           goto 70
 80       continue
       endif

!    ........................ computing kinetic energy matrix elements
!    in m-scheme
     !initialization
        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            Esp(jj,ii)=zero
         end do
        end do
!          write(110,*) '.... tk,ll,kk,Esp(ll,kk)....'

         write(*,'(a)',advance='no') &
     &  '   Transform M1B from J-scheme to M-scheme ...' 
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(HO,tnljm,ee,Esp,Ascale)
!!$OMP DO SCHEDULE(DYNAMIC)
        DO kk=1,HO%NLEV
        DO ll=1,HO%NLEV

          tk = tnljm%t(kk)  ! -1(p); +1 (n) 
          nk = tnljm%n(kk)  ! 1,2,3,...
          twomk = tnljm%twom(kk)
          ljk   = tnljm%lj(kk)

          tl = tnljm%t(ll)  ! -1(p); +1 (n) 
          nl = tnljm%n(ll)  ! 1,2,3,...
          twoml = tnljm%twom(ll)
          ljl   = tnljm%lj(ll)
          if(tk .ne. tl .or. ljk.ne.ljl .or. twomk .ne. twoml) cycle
          Esp(ll,kk) = ee(tk,nl,nk,ljk)*Ascale
        END DO
        END DO
!!$OMP END DO
!!$OMP END PARALLEL     
        write(*,*) ' done'
        end subroutine
