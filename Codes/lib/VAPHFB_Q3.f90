!    .................................................
        subroutine q3mume(Q30t)
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        double precision rnla

        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q30t(HO%NLEV,HO%NLEV)


        b3=HO%b_osc**3
        
!initializing

        do ii=1,HO%NLEV
        do jj=1,HO%NLEV
         Q30t(jj,ii)=zero
        end do
        end do

!        tnljm%t(kk) = mtisos(kk)
!          tnljm%n(kk) = nlindex(kk)
!          tnljm%l(kk) = lang(kk)
!          tnljm%twoj(kk) = jang(kk)
!          tnljm%lj(kk)   = (2*lang(kk)+jang(kk)-1)/2
!          tnljm%twom(kk) = mjang(kk)

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

	  delta_isos=0.d0
	  if (jta.eq.jtb) delta_isos=1.d0

	  factor=(((2.d0*jlb+1.d0)*7.d0)/(4.d0*pi*(2.d0*jla+1.d0)))**.5

	  call CJJ(6,2*jlb,2*jla,0,0,0,cbcomm)
	  
	  radial=b3*rnla(jnla,jla,3,jnlb,jlb) !Rnl(jnla,jnlb)

          suma_m2=0.d0
          suma_m1=0.d0
          suma_0 =0.d0
          suma_p1=0.d0
          suma_p2=0.d0

          DO MLA=0,2*jla
	   MMLA=-jla+MLA
           DO MLB=0,2*jlb
           MMLB=-jlb+MLB
	   
            DO IS=0,2,2   ! spin part
	     jms=-1+IS
	      
             call CJJ(2*jla,1,jja,2*MMLA,jms,jma,cb1)
             call CJJ(2*jlb,1,jjb,2*MMLB,jms,jmb,cb2)
	   

	     call CJJ(6,2*jlb,2*jla, 0,2*mmlb,2*mmla,cb_0)
 
 	     suma_0 =suma_0 +cb1*cb2*cb_0
            END DO
           END DO    
          END DO

!        jmy: the factor sqrt(7./(4.*pi)) is divided in which case
!             beta_L = 4*pi/(3AR^L) * Q_L
           Q30t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 !/sqrt(7./(4.*pi))

          if(abs(Q30t(ia,ib)).gt.1.d-5) write(193,'(2i4,f10.5)') ia,ib,Q30t(ia,ib)
         end do
        end do
        end subroutine

!    ......................................................
        subroutine q3mume_KB(Q20t)
!    ......................................................
!       Baranger and Kumar defintion:
!       M. Baranger and K. Kumar, Nucl. Phys. A 110, 490 (1968). 
!       Q/b^2 = sqrt(16*pi/5)[e_n*<r^2*Y_2> /b^2 + e_p *<r^2 Y_2>/b^2],
!       where, e_n =2N/A, e_p = 2Z/A, to be checked ... 
!
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

        double precision rnla

        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q30t(HO%NLEV,HO%NLEV)

        b3   = 1.d0
!       ...................................... red_fac = (3+1.5)/(3+2.5)
!        fac  = 2.d0/(sqrt(16*pi/5.d0))
        fac  = sqrt(7.d0/(4*pi))

!       ................................................
!        write(*,*) nucleon(0),nucleon(1),nucleon(2),facn,facp
!initializing

        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
         Q30t(jj,ii)=zero
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

          delta_isos=0.d0
          if (jta.eq.jtb )  delta_isos = fac       ! iso 

          factor=(((2.d0*jlb+1.d0)*7.d0)/(4.d0*pi*(2.d0*jla+1.d0)))**.5

          call CJJ(6,2*jlb,2*jla,0,0,0,cbcomm)

          radial=rnla(jnla,jla,2,jnlb,jlb)*b3 !Rnl(jnla,jnlb)

          suma_0 =0.d0

          DO MLA=0,2*jla
           MMLA=-jla+MLA
           DO MLB=0,2*jlb
           MMLB=-jlb+MLB

            DO IS=0,2,2
             jms=-1+IS

             call CJJ(2*jla,1,jja,2*MMLA,jms,jma,cb1)
             call CJJ(2*jlb,1,jjb,2*MMLB,jms,jmb,cb2)


             call CJJ(6,2*jlb,2*jla, 0,2*mmlb,2*mmla,cb_0)

             suma_0 =suma_0 +cb1*cb2*cb_0
            END DO
           END DO
          END DO

!             beta_L = 4*pi/(3AR^L) * Q_L
           Q30t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 /sqrt(7./(4.*pi))


          if(abs(Q30t(ia,ib)).gt.1.d-5) write(193,'(2i4,f10.5)') ia,ib,Q30t(ia,ib)

         end do
        end do

        end subroutine
