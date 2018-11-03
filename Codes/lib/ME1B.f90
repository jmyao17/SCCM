!    .................................................
        subroutine Rsq_ME(r2)
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)
        double precision rnla

        DIMENSION r2(HO%NLEV,HO%NLEV)
!initializing
        b2=HO%b_osc**2
        do ia=1,HO%NLEV
        do ib=1,HO%NLEV
           r2(ib,ia)=zero
        end do
        end do

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

         if (jja.eq.jjb .and. jma.eq.jmb .and. jta.eq.jtb) then 
           r2(ia,ib) = b2*rnla(jnla,jla,2,jnlb,jlb) 
         endif

         end do
        end do
        end subroutine
