        subroutine Reduced_ME1B()
        USE VAPHFB_PAR
        implicit none
        integer Lam,jja,jla,jta,jjb,jlb,jtb
        real*8 delta_isos,fact,b2,radial,reduce
        double precision rnla
        REAL(DP) :: Wigner_3j
        integer jna,jnb,lja,ljb,tnlja,tnljb       
        integer    SPB_twoj_lj,SPB_l_lj
 
        if(.NOT. ALLOCATED(rME1B%T)) allocate(rME1B%T(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2))

        do Lam=0,2,2

           b2=HO%b_osc**Lam
        do tnlja=1,FSPB%tnlj_fmax
           jta = FSPB%t(tnlja) 
           jna = FSPB%n(tnlja) ! starting from 1 
           lja = FSPB%lj(tnlja)
           jla = SPB_l_lj(lja)       ! NOT doubled 
           jja = SPB_twoj_lj(lja)    ! doubled
 
        do tnljb=1,FSPB%tnlj_fmax
           jtb = FSPB%t(tnljb) ! n(0), p(1)
           jnb = FSPB%n(tnljb) ! starting from 1 
           ljb = FSPB%lj(tnljb)
           jlb = SPB_l_lj(ljb)     ! NOT doubled 
           jjb = SPB_twoj_lj(ljb)    ! doubled
      
          if (jta.ne.jtb ) cycle
          delta_isos=0.d0
          ! only proton 
          if (jta.eq.1)   delta_isos=-1.d0
          fact = delta_isos*(iv(jla+jlb+Lam)+1.d0)/2.d0 &
      &          *sq(jja+1)/sq(2*Lam+1)

          radial=b2*rnla(jna,jla,Lam,jnb,jlb) !Rnl(jnla,jnlb)
          rME1B%T(tnlja,tnljb,Lam)=radial*sqrt((2*Lam+1)/(4.*pi))&
      &                     *fact*iv((jjb+2*Lam-1)/2)              &
      &                     *sqrt((jja+1.)*(jjb+1.))          &
      &                     *Wigner_3j(jja,jjb,2*Lam,1,-1,0)

       !   if(abs(rME1B%T(tnlja,tnljb,Lam)).gt. CHOP) &
       !   write(29,'(3i6,f15.8)') tnlja,tnljb,Lam,rME1B%T(tnlja,tnljb,Lam)
        end do
        end do
        enddo ! Lam
        return

        end
!      .................................................
        SUBROUTINE ZQ2(ZQME,ZRO,zQ_0p,zQ_0n,NDIM)
!      .................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)

        ZAUX(NDIM,NDIM) = zzero

!        do ii=1,NDIM
!        do jj=1,NDIM
!           ZAUX(ii,ii) =  ZAUX(ii,ii) + ZQME(ii,jj)*ZRO(jj,ii)
!        enddo
!        enddo
       
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)

        zQ_0p=zzero
        zQ_0n=zzero
        do ii=1,NDIM/2
         zQ_0p=zQ_0p+ZAUX(ii,ii)
         zQ_0n=zQ_0n+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE

!      .................................................
        SUBROUTINE ZQ1B(ZQME,ZRO,zQp,zQn,NDIM)
!      .................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)

        ZAUX(NDIM,NDIM) = zzero

!        do ii=1,NDIM
!        do jj=1,NDIM
!           ZAUX(ii,ii) =  ZAUX(ii,ii) + ZQME(ii,jj)*ZRO(jj,ii)
!        enddo
!        enddo

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)

        zQp=zzero
        zQn=zzero
        do ii=1,NDIM/2
         zQp=zQp+ZAUX(ii,ii)
         zQn=zQn+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE

!      .................................................
        SUBROUTINE ZQ2t(ZQME,ZRO,zQ_0p,zQ_0n,NDIM)
!      .................................................
!       transponent of Q2
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)


        call ZGEMM ('c','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)
        zQ_0p=zzero
        zQ_0n=zzero
        do ii=1,NDIM/2
         zQ_0p=zQ_0p+ZAUX(ii,ii)
         zQ_0n=zQ_0n+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE





        subroutine SingleParticleME()
        use vaphfb_par
        implicit none
        integer i,j
!        DIMENSION akin(NLEV,NLEV), zakin(NLEV,NLEV)
!       DIMENSION ZAZme(NLEV**2),ZANme(NLEV**2)
!        DIMENSION AZme(NLEV**2),ANme(NLEV**2)
!        DIMENSION AJZ_ME(NLEV2),AJZ2_ME(NLEV2)

!   ......................... allocate memory

      if(.NOT. ALLOCATED(cME1B%Q2_2t))  ALLOCATE(cME1B%Q2_2t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q2_1t))  ALLOCATE(cME1B%Q2_1t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q20t))  ALLOCATE(cME1B%Q20t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q21t))  ALLOCATE(cME1B%Q21t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q22t))  ALLOCATE(cME1B%Q22t(1:HO%NLEV*HO%NLEV))

      cME1B%Q2_2t = zero
      cME1B%Q2_1t = zero
      cME1B%Q20t = zero
      cME1B%Q21t = zero
      cME1B%Q22t = zero

!      if(.NOT. ALLOCATED(cME1B%Q2_2m))  ALLOCATE(cME1B%Q2_2m(1:HO%NLEV*HO%NLEV))
!      if(.NOT. ALLOCATED(cME1B%Q2_1m))  ALLOCATE(cME1B%Q2_1m(1:HO%NLEV*HO%NLEV))
!      if(.NOT. ALLOCATED(cME1B%Q20m))  ALLOCATE(cME1B%Q20m(1:HO%NLEV*HO%NLEV))
!      if(.NOT. ALLOCATED(cME1B%Q21m))  ALLOCATE(cME1B%Q21m(1:HO%NLEV*HO%NLEV))
!      if(.NOT. ALLOCATED(cME1B%Q22m))  ALLOCATE(cME1B%Q22m(1:HO%NLEV*HO%NLEV))


      if(.NOT. ALLOCATED(cME1B%zJxME))  ALLOCATE(cME1B%zJxME(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zJyME))  ALLOCATE(cME1B%zJyME(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zJzME))  ALLOCATE(cME1B%zJzME(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zJx2ME))  ALLOCATE(cME1B%zJx2ME(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zJy2ME))  ALLOCATE(cME1B%zJy2ME(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zJz2ME))  ALLOCATE(cME1B%zJz2ME(1:HO%NLEV*HO%NLEV))

      if(.NOT. ALLOCATED(cME1B%r2)) ALLOCATE(cME1B%r2(1:HO%NLEV*HO%NLEV))
!  ........... radii
       call Rsq_ME(cME1B%r2)

      if(Const%iQBType.eq.0) then
!  ................. Kumar-Baranger definition for Qn+Qp
        call q2mume_KB(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)

       

!  ................. standard defintion for Qn+Qp
       else if(Const%iQBType.eq.1) then
          call q2mume(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)

!         call q2mume_ME1B(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)
!  ................. read matrix elements from Nathan
       else if(Const%iQBType.eq.2) then
!         call Generate_ME1B4Q2(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)
       else
          stop 'E2 matrix elements are not set properly !'
       endif
!  ................. standard definition for Qn-Qp 
!          call q2mume_IsoV(cME1B%Q2_2m,cME1B%Q2_1m,cME1B%Q20m,cME1B%Q21m,cME1B%Q22m)

!       ..........angular momentum operators
       if(.NOT. ALLOCATED(H%zME1BM)) ALLOCATE(H%zME1BM(1:HO%NLEV,1:HO%NLEV))
        call Real2Complex(H%ME1BM,H%zME1BM,HO%NLEV)
!particle numbers
        CALL Complex_ZNme()
!angular momentum
!        write(*,*) ' computing s.p. elements for Jx,Jy,Jz'
        CALL zJiME(cME1B%ZJxME,cME1B%ZJyME,cME1B%zJzME,cME1B%ZJx2ME,cME1B%ZJy2ME,cME1B%zJz2ME,HO%NLEV)

!        do i=1,HO%NLEV
!        do j=1,HO%NLEV
!           if(abs(cME1B%ZJx2ME(i+j*(HO%NLEV-1))).gt.1.d-4) &
!      &     write(700,'(2i3,f10.5)') i,j,dreal(cME1B%ZJx2ME(i+j*(HO%NLEV-1)))
!        enddo
!        enddo

        return
        end

      SUBROUTINE ZN_0(ZROBP,Zprot_0,Zneut_0,N)
      implicit real*8(a-h,o-y)
      implicit complex*16 (z)

      Parameter (Zone=(1.d0,0.d0))
      Parameter (Zzero=(0.d0,0.d0))
      Parameter (Zimag=(0.d0,1.d0))

      Dimension ZROBP(N,N)

      zprot_0=zzero
      zneut_0=zzero

      do ii=1,N/2
       zprot_0=zprot_0+ZROBP(ii,ii)
       zneut_0=zneut_0+ZROBP(ii+N/2,ii+N/2)
      end do


      end subroutine



!       ...................................
        subroutine Parity_ip()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

      if(.NOT. ALLOCATED(ipme))  ALLOCATE(ipme(1:HO%NLEV,1:HO%NLEV))

        ipme = 0
        do ia=1,HO%NLEV
          do ib=1,HO%NLEV
             ila =tnljm%l(ia) ! lang(ia)
             ilb =tnljm%l(ib) ! lang(ia)
             if(iv(ila+ilb).eq.0) then 
                ipme(ib,ia)=iv(ila)
             else
                ipme(ib,ia)= 0
             endif 
          end do
       end do
       return
       end subroutine

!       ...................................
        subroutine Isospin_t3()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

      if(.NOT. ALLOCATED(t3me))  ALLOCATE(t3me(1:HO%NLEV,1:HO%NLEV))

        do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             t3me(jj,ii)=zero
             if((ii.eq.jj).AND.(ii.le.(HO%NLEV/2))) then
                t3me(jj,ii)=-1.d0
             else if((ii.eq.jj).AND.(ii.gt.(HO%NLEV/2))) then
                t3me(jj,ii)=1.d0
             end if
          end do
       end do
       end subroutine

!       ...................................
        subroutine ZNme()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
 
      if(.NOT. ALLOCATED(cME1B%AZme))  ALLOCATE(cME1B%AZme(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%ANme))  ALLOCATE(cME1B%ANme(1:HO%NLEV*HO%NLEV))

       ll=0
        do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             ll=ll+1
             cME1B%AZme(ll)=zero
             cME1B%ANme(ll)=zero
             if((ii.eq.jj).AND.(ii.le.(HO%NLEV/2))) then
                cME1B%AZme(ll)=one
             else if((ii.eq.jj).AND.(ii.gt.(HO%NLEV/2))) then
                cME1B%ANme(ll)=one
             end if
          end do
       end do
       end subroutine

!       ...................................
        subroutine Complex_ZNme()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

      if(.NOT. ALLOCATED(cME1B%zAZme))  ALLOCATE(cME1B%zAZme(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%zANme))  ALLOCATE(cME1B%zANme(1:HO%NLEV*HO%NLEV))

       ll=0
        do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             ll=ll+1
             cME1B%zAZme(ll)=zzero
             cME1B%zANme(ll)=zzero
             if((ii.eq.jj).AND.(ii.le.(HO%NLEV/2))) then
                cME1B%zAZme(ll)=zone
             else if((ii.eq.jj).AND.(ii.gt.(HO%NLEV/2))) then
                cME1B%zANme(ll)=zone
             end if
          end do
       end do
       end subroutine

        subroutine zJiME(ZAJX,ZAJY,ZAJZ,ZAJX2,ZAJY2,ZAJZ2,NLEV)
        USE VAPHFB_Par
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZAJX(NLEV,NLEV)
        DIMENSION ZAJX2(NLEV,NLEV)
        DIMENSION ZAJY(NLEV,NLEV)
        DIMENSION ZAJY2(NLEV,NLEV)
        DIMENSION ZAJZ(NLEV,NLEV)
        DIMENSION ZAJZ2(NLEV,NLEV)

        do ii=1,HO%NLEV
         do jj=1,HO%NLEV
         ZAJX(jj,ii)=zzero
         ZAJX2(jj,ii)=zzero
         ZAJY(jj,ii)=zzero
         ZAJY2(jj,ii)=zzero
         ZAJZ(jj,ii)=zzero
         ZAJZ2(jj,ii)=zzero
         end do
        end do

        do ia=1,HO%NLEV
         ina=tnljm%n(ia) ! nlindex(ia)
         ila =tnljm%l(ia) ! lang(ia)
         ija =tnljm%twoj(ia) ! jang(ia)
         ima =tnljm%twom(ia) ! mjang(ia)
         ita =tnljm%t(ia)    ! mtisos(ia)
         rja=ija/2.d0
         rma=ima/2.d0

         do ib=1,HO%NLEV
         inb=tnljm%n(ib) ! nlindex(ib)
         ilb =tnljm%l(ib) ! lang(ib)
         ijb =tnljm%twoj(ib) ! jang(ib)
         imb =tnljm%twom(ib) ! mjang(ib)
         itb =tnljm%t(ib)    ! mtisos(ib)
         rjb=ijb/2.d0
         rmb=imb/2.d0


         fp1=rjb*(rjb+1.d0)-rmb*(rmb+1.d0)  ! j(j+1) - m(m+1)
         fm1=rjb*(rjb+1.d0)-rmb*(rmb-1.d0)  ! j(j+1) - m(m-1)

         fp1p2=rjb*(rjb+1.d0)-(rmb+1.d0)*(rmb+2.d0) ! j(j+1)-(m+1)(m+2)
         fm1m2=rjb*(rjb+1.d0)-(rmb-1.d0)*(rmb-2.d0) ! j(j+1)-(m-1)(m-2)
         f02=  rjb*(rjb+1.d0)-rmb**2                ! j(j+1)-2m


         if ((ina.eq.inb).AND.(ila.eq.ilb)           &
     &   .AND.(ija.eq.ijb).AND.(ita.eq.itb)) then


         delta_ma_mbp1=0.d0
         delta_ma_mbm1=0.d0
         delta_ma_mb  =0.d0
         delta_ma_mbp2=0.d0
         delta_ma_mbm2=0.d0

         if ((ima.eq.(imb+2))) delta_ma_mbp1=1.d0
         if ((ima.eq.(imb-2))) delta_ma_mbm1=1.d0
         if ((ima.eq.(imb+4))) delta_ma_mbp2=1.d0
         if ((ima.eq.(imb-4))) delta_ma_mbm2=1.d0
         if ((ima.eq.imb)) delta_ma_mb=1.d0

         ZAJX(ia,ib)=0.5d0*(dsqrt(fp1)*delta_ma_mbp1+     &
     &   dsqrt(fm1)*delta_ma_mbm1)*zone

         ZAJY(ia,ib)=-0.5d0*(dsqrt(fp1)*delta_ma_mbp1-    &
     &   dsqrt(fm1)*delta_ma_mbm1)*zimag

         ZAJZ(ia,ib)=rmb*delta_ma_mb*zone

         ZAJX2(ia,ib)=0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+  &
     &   dsqrt(fm1*fm1m2)*delta_ma_mbm2+                       &
     &   2*f02*delta_ma_mb)*zone

         ZAJY2(ia,ib)=-0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+  &
     &   dsqrt(fm1*fm1m2)*delta_ma_mbm2-                        &
     &   2*f02*delta_ma_mb)*zone

         ZAJZ2(ia,ib)=zone*(rmb*delta_ma_mb)**2

         end if

         end do
        end do


        end subroutine



       subroutine JZ_JZ2_ME(AJX,AJY,AJZ,AJX2,AJY2,AJZ2,NLEV)
     
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

        DIMENSION AJX(NLEV,NLEV)
        DIMENSION AJX2(NLEV,NLEV)
        DIMENSION AJY(NLEV,NLEV)
        DIMENSION AJY2(NLEV,NLEV)
        DIMENSION AJZ(NLEV,NLEV)
        DIMENSION AJZ2(NLEV,NLEV)
!        b2=HO%b_osc**2
!initializing
         do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            AJX(jj,ii)=zero
            AJX2(jj,ii)=zero
            AJY(jj,ii)=zero
            AJY2(jj,ii)=zero
            AJZ(jj,ii)=zero
            AJZ2(jj,ii)=zero
        end do
        end do

        do ia=1,NLEV
         ina=tnljm%n(ia) ! nlindex(ia)
         ila =tnljm%l(ia) ! lang(ia)
         ija =tnljm%twoj(ia) ! jang(ia)
         ima =tnljm%twom(ia) ! mjang(ia)
         ita =tnljm%t(ia)    ! mtisos(ia)

         rja=ija/2.d0
         rma=ima/2.d0

         do ib=1,NLEV
         inb=tnljm%n(ib) ! nlindex(ia)
         ilb =tnljm%l(ib) ! lang(ia)
         ijb =tnljm%twoj(ib) ! jang(ia)
         imb =tnljm%twom(ib) ! mjang(ia)
         itb =tnljm%t(ib)    ! mtisos(ia)

         rjb=ijb/2.d0
         rmb=imb/2.d0

 
         fp1=rjb*(rjb+1.d0)-rmb*(rmb+1.d0)
         fm1=rjb*(rjb+1.d0)-rmb*(rmb-1.d0)

         fp1p2=rjb*(rjb+1.d0)-(rmb+1.d0)*(rmb+2.d0)
         fm1m2=rjb*(rjb+1.d0)-(rmb-1.d0)*(rmb-2.d0)
         f02=  rjb*(rjb+1.d0)-rmb**2

         if ((ina.eq.inb).AND.(ila.eq.ilb)            &
     &   .AND.(ija.eq.ijb).AND.(ita.eq.itb)) then
         
         delta_ma_mbp1=0.d0
         delta_ma_mbm1=0.d0
         delta_ma_mb  =0.d0
         delta_ma_mbp2=0.d0
         delta_ma_mbm2=0.d0

         if ((ima.eq.(imb+2))) delta_ma_mbp1=1.d0
         if ((ima.eq.(imb-2))) delta_ma_mbm1=1.d0
         if ((ima.eq.(imb+4))) delta_ma_mbp2=1.d0
         if ((ima.eq.(imb-4))) delta_ma_mbm2=1.d0
         if ((ima.eq.imb)) delta_ma_mb=1.d0
 
         AJX(ia,ib)=0.5d0*(dsqrt(fp1)*delta_ma_mbp1+               &
     &                    dsqrt(fm1)*delta_ma_mbm1)  
         AJY(ia,ib)=0.5d0*(dsqrt(fp1)*delta_ma_mbp1-               &
     &                     dsqrt(fm1)*delta_ma_mbm1)
     
         AJZ(ia,ib)=rmb*delta_ma_mb
     
         AJX2(ia,ib)=0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+       &
     &                       dsqrt(fm1*fm1m2)*delta_ma_mbm2+         &
     &                                2*f02*delta_ma_mb)

         AJY2(ia,ib)=0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+       &
     &                      dsqrt(fm1*fm1m2)*delta_ma_mbm2-       &
     &                                2*f02*delta_ma_mb)
         
             AJZ2(ia,ib)=(rmb*delta_ma_mb)**2
           end if
          end do
        end do
        end subroutine



        SUBROUTINE Qpair_0(QME,akapa10,akapa01,Q_0,NDIM)
        
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
!        Parameter (zero=0.d0)
!        Parameter (one=1.d0)
        
        DIMENSION QME(NDIM,NDIM)
        DIMENSION akapa10(NDIM,NDIM),akapa01(NDIM,NDIM)
        DIMENSION AUX1(NDIM,NDIM)
        DIMENSION AUX2(NDIM,NDIM)

        
        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,    &
     &  QME,NDIM,akapa10,NDIM,zero,AUX1,NDIM)   

        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,    &
     &  QME,NDIM,akapa01,NDIM,zero,AUX2,NDIM)   
        
        
        Q_0=zero
        do ii=1,NDIM
         Q_0=Q_0+(AUX1(ii,ii)+AUX2(ii,ii))/2.d0
        end do
        
        END SUBROUTINE
