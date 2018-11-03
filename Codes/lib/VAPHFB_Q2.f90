
!    .................................................
        subroutine q2mume_ME1B(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
!    .................................................
!     Book by J. Suhonen Eq.(6.6) and (6.23) with condon-shortley phase
!    .................................................

        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        double precision rnla
        REAL(DP) :: CG
        REAL(DP) :: Wigner_3j 


        DIMENSION Rnl(HO%nmax,HO%nmax) ! (2,2)
        DIMENSION Q2_2t(HO%NLEV,HO%NLEV)
        DIMENSION Q2_1t(HO%NLEV,HO%NLEV)
        DIMENSION Q20t(HO%NLEV,HO%NLEV)
        DIMENSION Q21t(HO%NLEV,HO%NLEV)
        DIMENSION Q22t(HO%NLEV,HO%NLEV)

!        open(193,file='Suhonen_Q20_'//Input%chwHO//'_'// &
!      &     Input%cValID//'_me1b.dat',status='unknown')

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

        open(191,file='Q20_'//Input%chwHO//'_'// &
      &     Input%cValID//'_me1b.dat',status='unknown')


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

!             beta_L = 4*pi/(3AR^L) * Q_L
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2i !/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1  !/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0  !/sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1 !/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2 !/sqrt(5./(4.*pi))

!          if(jta.eq.-1 .and. abs(Q20t(ia,ib)).gt.1.d-5) &
!      &   write(191,'(2i4,f10.5)') ia,ib,Q20t(ia,ib)
         end do
        end do

        return
!        DO ll=1,HO%NLEV
!        DO kk=1,HO%NLEV
!          if(tnljm%t(ll).eq.-1 .and. abs(Q20t(ll,kk)).gt.1.d-5) &
!       &  write(191,'(2i4,f10.5)') ll,kk,Q20t(ll,kk) !f1bJ(tk,nl,nk,ljl,ljk)
!        enddo
!        enddo

        end subroutine

!    ......................................................
        subroutine q2mume_IsoV(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
!    ......................................................
!       <a | Q_2mu | b>/ sqrt(5/4pi)
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
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2 !/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1 !/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0  !/sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1 !/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2 !/sqrt(5./(4.*pi))

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

        b2   = 1.d0
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

           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0 /sqrt(5./(4.*pi))
           Q21t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p1/sqrt(5./(4.*pi))
           Q22t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_p2/sqrt(5./(4.*pi))


          if(abs(Q20t(ia,ib)).gt.1.d-5) write(191,'(2i4,f10.5)') ia,ib,Q20t(ia,ib)

         end do
        end do

        end subroutine



!    ......................................................
        subroutine q2mumeIso_KB(Q2_2t,Q2_1t,Q20t,Q21t,Q22t)
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

        b2   = 1.d0
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
          if (jta.eq.jtb .and. jta.eq.1)  delta_isos= fac       !isospin should be the same
          if (jta.eq.jtb .and. jta.eq.-1) delta_isos=-fac       !


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

!        the factor sqrt(5./(4.*pi)) is divided in which case
           Q2_2t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m2/sqrt(5./(4.*pi))
           Q2_1t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_m1/sqrt(5./(4.*pi))
           Q20t(ia,ib) =delta_isos*factor*cbcomm*radial*suma_0/sqrt(5./(4.*pi))
           Q21t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_p1/sqrt(5./(4.*pi))
           Q22t(ia,ib)=delta_isos*factor*cbcomm*radial*suma_p2/sqrt(5./(4.*pi))


!          if(abs(Q20t(ia,ib)).gt.1.d-5) write(191,'(2i4,f10.5)') ia,ib,Q20t(ia,ib)

         end do
        end do

        end subroutine
