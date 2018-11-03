!    .................................................
        subroutine q2mume(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        double precision rnla

        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q2_2t(HO%NLEV,HO%NLEV)
        DIMENSION Q2_1t(HO%NLEV,HO%NLEV)
        DIMENSION Q20t(HO%NLEV,HO%NLEV)
        DIMENSION Q21t(HO%NLEV,HO%NLEV)
        DIMENSION Q22t(HO%NLEV,HO%NLEV)


        b2=HO%b_osc**2
        
!initializing

        do ii=1,HO%NLEV
        do jj=1,HO%NLEV
         Q2_2t(jj,ii)=zero
         Q2_1t(jj,ii)=zero
         Q20t(jj,ii)=zero
         Q21t(jj,ii)=zero
         Q22t(jj,ii)=zero
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

	  factor=(((2.d0*jlb+1.d0)*5.d0)/(4.d0*pi*(2.d0*jla+1.d0)))**.5

	  call CJJ(4,2*jlb,2*jla,0,0,0,cbcomm)
	  
	  radial=b2*rnla(jnla,jla,2,jnlb,jlb) !Rnl(jnla,jnlb)

          suma_m2=0.d0
          suma_m1=0.d0
          suma_0 =0.d0
          suma_p1=0.d0
          suma_p2=0.d0

          DO MLA=0,2*jla
	   MMLA=-jla+MLA
           DO MLB=0,2*jlb
           MMLB=-jlb+MLB
	   
            DO IS=0,2,2
	     jms=-1+IS
	      
             call CJJ(2*jla,1,jja,2*MMLA,jms,jma,cb1)
             call CJJ(2*jlb,1,jjb,2*MMLB,jms,jmb,cb2)
	   

	     call CJJ(4,2*jlb,2*jla,-4,2*mmlb,2*mmla,cbm2)
	     call CJJ(4,2*jlb,2*jla,-2,2*mmlb,2*mmla,cbm1)
	     call CJJ(4,2*jlb,2*jla, 0,2*mmlb,2*mmla,cb_0)
	     call CJJ(4,2*jlb,2*jla, 2,2*mmlb,2*mmla,cbp1)
	     call CJJ(4,2*jlb,2*jla, 4,2*mmlb,2*mmla,cbp2)
 
 	     suma_m2=suma_m2+cb1*cb2*cbm2
 	     suma_m1=suma_m1+cb1*cb2*cbm1
 	     suma_0 =suma_0 +cb1*cb2*cb_0
 	     suma_p1=suma_p1+cb1*cb2*cbp1
 	     suma_p2=suma_p2+cb1*cb2*cbp2
            END DO
           END DO    
          END DO

!        jmy: the factor sqrt(5./(4.*pi)) is divided in which case
!             beta_L = 4*pi/(3AR^L) * Q_L
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 /sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2/sqrt(5./(4.*pi))

         end do
        end do
        end subroutine

!    ......................................................
        subroutine q2mume_IsoV(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
!    ......................................................

        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

        double precision rnla

        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q2_2t(HO%NLEV,HO%NLEV)
        DIMENSION Q2_1t(HO%NLEV,HO%NLEV)
        DIMENSION Q20t(HO%NLEV,HO%NLEV)
        DIMENSION Q21t(HO%NLEV,HO%NLEV)
        DIMENSION Q22t(HO%NLEV,HO%NLEV)

        b2=HO%b_osc**2

!initializing

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
          if (jta.eq.jtb .and. jta.eq.1)  delta_isos= 1.d0       ! isospin should be the same
          if (jta.eq.jtb .and. jta.eq.-1) delta_isos=-1.d0       ! isospin should be the same

          factor=(((2.d0*jlb+1.d0)*5.d0)/(4.d0*pi*(2.d0*jla+1.d0)))**.5

          call CJJ(4,2*jlb,2*jla,0,0,0,cbcomm)
        
          radial=b2*rnla(jnla,jla,2,jnlb,jlb) !Rnl(jnla,jnlb)

          suma_m2=0.d0
          suma_m1=0.d0
          suma_0 =0.d0
          suma_p1=0.d0
          suma_p2=0.d0

          DO MLA=0,2*jla
           MMLA=-jla+MLA
           DO MLB=0,2*jlb
           MMLB=-jlb+MLB
        
            DO IS=0,2,2
             jms=-1+IS
        
             call CJJ(2*jla,1,jja,2*MMLA,jms,jma,cb1)
             call CJJ(2*jlb,1,jjb,2*MMLB,jms,jmb,cb2)
        

             call CJJ(4,2*jlb,2*jla,-4,2*mmlb,2*mmla,cbm2)
             call CJJ(4,2*jlb,2*jla,-2,2*mmlb,2*mmla,cbm1)
             call CJJ(4,2*jlb,2*jla, 0,2*mmlb,2*mmla,cb_0)
             call CJJ(4,2*jlb,2*jla, 2,2*mmlb,2*mmla,cbp1)
             call CJJ(4,2*jlb,2*jla, 4,2*mmlb,2*mmla,cbp2)

             suma_m2=suma_m2+cb1*cb2*cbm2
             suma_m1=suma_m1+cb1*cb2*cbm1
             suma_0 =suma_0 +cb1*cb2*cb_0
             suma_p1=suma_p1+cb1*cb2*cbp1
             suma_p2=suma_p2+cb1*cb2*cbp2
            END DO
           END DO
          END DO

!        jmy: the factor sqrt(5./(4.*pi)) is divided in which case
!             beta_L = 4*pi/(3AR^L) * Q_L
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 /sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2/sqrt(5./(4.*pi))

         end do
        end do
        end subroutine


!    ......................................................
        subroutine q2mume_KB(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
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
        DIMENSION Q2_2t(HO%NLEV,HO%NLEV)
        DIMENSION Q2_1t(HO%NLEV,HO%NLEV)
        DIMENSION Q20t(HO%NLEV,HO%NLEV)
        DIMENSION Q21t(HO%NLEV,HO%NLEV)
        DIMENSION Q22t(HO%NLEV,HO%NLEV)

        b2=HO%b_osc**2
        b2   = 4.2627833 !1.d0
!       ...................................... red_fac = (3+1.5)/(3+2.5)
!        fac  = 2.d0/(sqrt(16*pi/5.d0))
        fac  = sqrt(5.d0/(4*pi))

!       ................................................
!        write(*,*) nucleon(0),nucleon(1),nucleon(2),facn,facp
!initializing

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

          factor=(((2.d0*jlb+1.d0)*5.d0)/(4.d0*pi*(2.d0*jla+1.d0)))**.5

          call CJJ(4,2*jlb,2*jla,0,0,0,cbcomm)

          radial=rnla(jnla,jla,2,jnlb,jlb)*b2 !Rnl(jnla,jnlb)

          suma_m2=0.d0
          suma_m1=0.d0
          suma_0 =0.d0
          suma_p1=0.d0
          suma_p2=0.d0

          DO MLA=0,2*jla
           MMLA=-jla+MLA
           DO MLB=0,2*jlb
           MMLB=-jlb+MLB

            DO IS=0,2,2
             jms=-1+IS

             call CJJ(2*jla,1,jja,2*MMLA,jms,jma,cb1)
             call CJJ(2*jlb,1,jjb,2*MMLB,jms,jmb,cb2)


             call CJJ(4,2*jlb,2*jla,-4,2*mmlb,2*mmla,cbm2)
             call CJJ(4,2*jlb,2*jla,-2,2*mmlb,2*mmla,cbm1)
             call CJJ(4,2*jlb,2*jla, 0,2*mmlb,2*mmla,cb_0)
             call CJJ(4,2*jlb,2*jla, 2,2*mmlb,2*mmla,cbp1)
             call CJJ(4,2*jlb,2*jla, 4,2*mmlb,2*mmla,cbp2)

             suma_m2=suma_m2+cb1*cb2*cbm2
             suma_m1=suma_m1+cb1*cb2*cbm1
             suma_0 =suma_0 +cb1*cb2*cb_0
             suma_p1=suma_p1+cb1*cb2*cbp1
             suma_p2=suma_p2+cb1*cb2*cbp2
            END DO
           END DO
          END DO

!        jmy: the factor sqrt(5./(4.*pi)) is divided in which case
!             beta_L = 4*pi/(3AR^L) * Q_L
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 /sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2/sqrt(5./(4.*pi))


          if(abs(Q20t(ia,ib)).gt.1.d-5) write(191,'(2i4,f10.5)') ia,ib,Q20t(ia,ib)

         end do
        end do

        end subroutine
