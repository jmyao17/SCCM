        ! from JT to M scheme 

         subroutine INT_ME2B_JT2M()

         USE VAPHFB_PAR
         use omp_lib
         implicit none

        integer itype
        integer ii,kk,mj,kk_n
        integer ia,ib,ic,id,na,lja,nb,ljb,nc,ljc,nd,ljd
        integer mjinc,jttot(1:2)
        integer jmax,nlj,iabcd
        real*8  hja,hjb,hjc,hjd
        integer it12,MMT,iJ,iT
        real*8  ajtot,attot,ame_coup
        real*8  phasJT,phasab,phascd,delta_ab,delta_cd
        integer jja,jma,jta,jna,jlja,jacoup
        integer jjb,jmb,jtb,jnb,jljb,jbcoup
        integer jjc,jmc,jtc,jnc,jljc,jccoup
        integer jjd,jmd,jtd,jnd,jljd,jdcoup
        integer JJ,MMJ,MT
        integer icheck
        real*8  AN_INV,AN_AB,AN_CD,cb1,cb2,cb3,cb4,suma,sumando
        real*8, dimension(:,:,:,:,:,:,:), allocatable :: V_AB_CD_JTMT


       open(1,file=trim(INT_DIR)//trim(File%f2bJTMT),status='old')
       open(99,file=trim(INT_DIR)//trim(File%f2bm),&
      &        form='unformatted',status='unknown')

       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))

6      read(1,*,end=7) ia,ib,ic,id,iJ,iT,MMT,ame_coup
       hja=SPB%twoj(ia)/2.d0
       hjb=SPB%twoj(ib)/2.d0
       hjc=SPB%twoj(ic)/2.d0
       hjd=SPB%twoj(id)/2.d0
       ajtot=1.d0*iJ ! jtot(iJ)
       attot=1.d0*iT ! jttot(iT)

       phasab=(-1.d0)**(hja+hjb+ajtot+attot)
       phascd=(-1.d0)**(hjc+hjd+ajtot+attot)

       V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,ic,id,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme
       V_AB_CD_JTMT(ib,ia,id,ic,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ia,ib,id,ic,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
!
       V_AB_CD_JTMT(ic,id,ia,ib,iJ,iT,MMT)= ame_coup*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ia,ib,iJ,iT,MMT)= ame_coup*phascd*H%ddd_tbme
       V_AB_CD_JTMT(id,ic,ib,ia,iJ,iT,MMT)= ame_coup*phasab*phascd*H%ddd_tbme
       V_AB_CD_JTMT(ic,id,ib,ia,iJ,iT,MMT)= ame_coup*phasab*H%ddd_tbme

       goto 6
  7   continue
      close(1)

!     ..............................................................................
!      OMP
!     ..............................................................................
      write(*,'(a30)') '     Transform NN from J to M: ' 

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(CG_Save,H,HO,tnljm,V_AB_CD_JTMT) 
!$OMP DO schedule(dynamic)  
      do iabcd = 1, H%iabcd_max
        ia = H%ka(iabcd)
        jja=tnljm%twoj(ia)  ! jang(ia)
        jma=tnljm%twom(ia)
        jta=tnljm%t(ia)      !mtisos(ia) 
        jna=tnljm%n(ia)
        jlja=tnljm%lj(ia)
        jacoup= tnljm%nlj(ia) ! jindex(ia)

        ib = H%kb(iabcd)
        jjb=tnljm%twoj(ib) ! jang(ib)
        jmb=tnljm%twom(ib)  ! mjang(ib)
        jtb=tnljm%t(ib)  !mtisos(ib)
        jnb=tnljm%n(ib)
        jljb=tnljm%lj(ib)
        jbcoup= tnljm%nlj(ib) ! jindex(ib)
        delta_ab=0.d0
        if(jna.eq.jnb .and. jlja.eq.jljb) delta_ab=1.d0
        ic = H%kc(iabcd)
         jjc=tnljm%twoj(ic)   ! jang(ic)
         jmc=tnljm%twom(ic)   ! mjang(ic)
         jtc=tnljm%t(ic)      ! mtisos(ic)
         jnc=tnljm%n(ic)
         jljc=tnljm%lj(ic)
         jccoup= tnljm%nlj(ic) !jindex(ic)
          id = H%kd(iabcd)
          jjd=tnljm%twoj(id)   !jang(id)
          jmd=tnljm%twom(id)   !mjang(id)
          jtd=tnljm%t(id)      !mtisos(id)
          jnd=tnljm%n(id)
          jljd=tnljm%lj(id)
          jdcoup= tnljm%nlj(id)  !jindex(id)
          delta_cd=0.d0
        if(jnc.eq.jnd .and. jljc.eq.jljd) delta_cd=1.d0

          suma=0.d0
          DO JJ=0,HO%twojmax
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ

              cb1  = CG_Save((jja+1)/2,(jma+1)/2,(jjb+1)/2,(jmb+1)/2,JJ,MMJ) ! J12,M12 are
              cb2  = CG_Save((jjc+1)/2,(jmc+1)/2,(jjd+1)/2,(jmd+1)/2,JJ,MMJ) ! J12,M12 are

           DO IT=0,HO%tmax

            phasJT=(-1.d0)**(JJ+IT)
!           ........................ normalization factor
            AN_AB=sqrt(1.d0-delta_ab*phasJT)/(1.d0+delta_ab)
            AN_CD=sqrt(1.d0-delta_cd*phasJT)/(1.d0+delta_cd)

            if(AN_AB.le.1e-15) cycle
            if(AN_CD.le.1e-15) cycle

         
            AN_INV=1.d0/(AN_AB*AN_CD)

            DO MT=0,2*IT
               MMT=-IT+MT

!              call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)
!              call CJJ(1,1,2*IT,jtc,jtd,2*MMT,cb4)

              cb3  = CG_Save(1,(jta+1)/2,1,(jtb+1)/2,IT,MMT) ! J12,M12 are NOT doubled 
              cb4  = CG_Save(1,(jtc+1)/2,1,(jtd+1)/2,IT,MMT) ! J12,M12 are NOT doubled 

!             if(abs(cg1-cb3).gt.CHOP) stop 'Error: in CG_Save '

             sumando=AN_INV*cb1*cb2*cb3*cb4                    &
     &        *V_AB_CD_JTMT(jacoup,jbcoup,jccoup,jdcoup,JJ,IT,MMT)
             suma=suma+sumando
            END DO
           END DO

          END DO
          END DO

           H%ME2BM(iabcd) = suma
       enddo ! iabcd         
!$OMP END DO 
!$OMP END PARALLEL
!     ..............................................................................
      write(*,*)
      write(*,'(a30)') ' Writing ME2B in m-scheme into file: '

       icheck = H%iabcd_max/100
       do iabcd = 1, H%iabcd_max
         if(mod(iabcd,icheck).eq.0) call progress_bar(int(iabcd/icheck))
         write(99) H%ME2BM(iabcd) ! suma
       enddo ! iabcd  
      write(*,*)
      
       return
        end
