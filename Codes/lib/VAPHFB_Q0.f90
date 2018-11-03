       subroutine q0_ME1B()
!    .................................................
!     Book by J. Suhonen Eq.(6.6) and (6.23) with condon-shortley phase
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        double precision rnla
        REAL(DP) :: CG
        REAL(DP) :: Wigner_3j 


        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)


        b2=HO%b_osc**2
        
!   ....... initializing

        do ii=1,HO%NLEV
        do jj=1,HO%NLEV
         Q2_2t(jj,ii)=zero
         Q2_1t(jj,ii)=zero
         Q20t(jj,ii)=zero
         Q21t(jj,ii)=zero
         Q22t(jj,ii)=zero
        end do
        end do

        do ia=1,HO%NLEV
         jnla=tnljm%n(ia) ! nlindex(ia)=1,2,..
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
          
  
          radial=b2*rnla(jnla,jla,2,jnlb,jlb) !Rnl(jnla,jnlb)
          reduce=radial*sqrt(5./(4.*pi))*iv((jjb-1)/2) &
      &                *sqrt((jja+1.)*(jjb+1.))          &
      &                *Wigner_3j(jja,jjb,2*2,1,-1,0)
         

           fac = iv(abs(jjb-jmb)/2)/sq(5)  ! 
           cg1=fac*CG(jja,jma,jjb,-jmb,2*2,-2*2)
           cg2=fac*CG(jja,jma,jjb,-jmb,2*2,-1*2)
           cg3=fac*CG(jja,jma,jjb,-jmb,2*2,0)
           cg4=fac*CG(jja,jma,jjb,-jmb,2*2,1*2)
           cg5=fac*CG(jja,jma,jjb,-jmb,2*2,2*2)
 
!        ...............................................................
!         beta_L = 4*pi/(3AR^L) * <r^L*Y_L0>
!        ...............................................................
           Q2_2t(ia,ib)=reduce*cg1 !/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=reduce*cg2 !/sqrt(5./(4.*pi))
           Q20t(ia,ib) =reduce*cg3 !/sqrt(5./(4.*pi))
           Q21t(ia,ib) =reduce*cg4 !/sqrt(5./(4.*pi))
           Q22t(ia,ib) =reduce*cg5 !/sqrt(5./(4.*pi))

!          if(jta.eq.-1 .and. abs(Q20t(ia,ib)).gt.1.d-5) &
!      &   write(193,'(2i4,2f10.5)') ia,ib,Q20t(ia,ib),reduce
         end do
        end do

        end subroutine
