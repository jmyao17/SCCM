!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine with the single particle energies expressed in the uncoupled basis n,l,m_l,m_s_mt!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine kinetic(Esp,Ascale)
        USE VAPHFB_PAR
        implicit none 
        Integer t1,n1,n2,LJ
        real*8 Ascale,ee(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax)
        Integer ii,jj,lk 
        Real*8 ee_t,Esp(HO%NLEV,HO%NLEV)
        integer kk,ll,tk,nk,twojk,twomk,ljk,tl,nl,twojl,twoml,ljl
        
!     ............. initialization
      ee(-1:1,HO%nmax,1:HO%nmax,0:HO%ljmax) = zero

       if(Input%IntIMSRG .ne. 1) then 
       write(*,*) 'Read S.P.E. from SPE_',Input%cIntID,'.dat' 
          open(21,file="../Int/SPE_"//Input%cIntID//".dat",status="old") 
       endif 

       if(Input%IntIMSRG .eq. 1) then
       write(*,*) 'Read S.P.E. from IMSRG_',Input%cFlow,'_SPE_',Input%cIntID,'.dat' 
       open(21,file="../Int/IMSRG_"//Input%cFlow//"_SPE_"//Input%cIntID//".dat",status="old")        
       endif
!     .................. read sp part of H
!     t1=0 (n); 1 (p); 2 (n+p)
       read(21,*) H%E0
 70    read(21,*,end=80) t1, LJ, n1, n2, ee_t
!          ee(iv(t1),n1+1,LJ) = ee_t        ! t1=0 (n); 1 (p)
         if(Input%IntType .eq. 0) ee(iv(t1),n1+1,n2+1,LJ) = ee_t              
         if(Input%IntType .eq. 1) ee(iv(t1),n1+1,n2+1,LJ) = ee_t*(one-one/Nucl%nucleon(2))
           write(*,*) t1, LJ, n1, n2, ee(iv(t1),n1+1,n2+1,LJ)
          goto 70
 80    continue

        !initialization
        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            Esp(jj,ii)=zero
         end do
        end do
          write(110,*) '.... tk,ll,kk,Esp(ll,kk)....'

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
          Esp(ll,kk) = ee(tk,nl,nk,ljk)*Ascale
          write(110,*) tk,ll,kk,Esp(ll,kk)
        END DO
        END DO
     
        end subroutine
