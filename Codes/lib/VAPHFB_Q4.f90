!    .................................................
        subroutine q40me_ME1B(Q40,iT)
!    .................................................
!      iT =0: isoscalar
!      iT =1: isovector
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        double precision rnla
        REAL(DP) :: CG
        REAL(DP) :: Wigner_3j 


        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q40(HO%NLEV,HO%NLEV)

        b4=HO%b_osc**4
        
!   ....... initializing

        do ii=1,HO%NLEV
        do jj=1,HO%NLEV
           Q40(jj,ii)=zero
        end do
        end do

        ! multipolarity
        Lambda =4
        Lambda2=2*Lambda+1
        do ia=1,HO%NLEV
         jnla=tnljm%n(ia) ! nlindex(ia)
         jla =tnljm%l(ia) ! lang(ia)
         jja =tnljm%twoj(ia) ! jang(ia)
         jma =tnljm%twom(ia) ! mjang(ia)
         jta =tnljm%t(ia)    ! mtisos(ia)
         do ib=1,HO%NLEV
         jnlb=tnljm%n(ib) ! nlindex(ib)
         jlb =tnljm%l(ib) ! lang(ib)
         jjb =tnljm%twoj(ib) ! jang(ib)
         jmb =tnljm%twom(ib) ! mjang(ib)
         jtb =tnljm%t(ib)    ! mtisos(ib)

          reduce = 0.d0
          if (jta.ne.jtb .or. iv(jla+jlb).ne.1) cycle 
          delta_isos= 1.d0
          if (iT.eq.1 .and. jta.eq.1)  delta_isos= 1.d0       !isospin should be the same
          if (iT.eq.1 .and. jta.eq.-1) delta_isos=-1.d0       !isospin should be the same          
  
          radial=b4*rnla(jnla,jla,Lambda,jnlb,jlb) !Rnl(jnla,jnlb)
          reduce=radial*sqrt(Lambda2/(4.*pi))*iv((jjb-1)/2) &
      &                *sqrt((jja+1.)*(jjb+1.))          &
      &                *Wigner_3j(jja,jjb,2*Lambda,1,-1,0)
         
           fac = iv(abs(jjb-jmb)/2)/sq(Lambda2)  ! 
           cg3 = fac*CG(jja,jma,jjb,-jmb,2*Lambda,0)
 
!        ...............................................................

           Q40(ia,ib) = delta_isos*reduce*cg3 
           write(300+iT,'(2i6,f12.8)') ia,ib,Q40(ia,ib)
         end do
        end do

        end subroutine
