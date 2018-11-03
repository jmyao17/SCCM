        subroutine QLreduced_ME1B()
!       .......................................
!       <Jf,Kf || Q_L || Ji,Ki>
!       where
!       - Q_2 = r^2 Y_2
!       - Q_0 = r^2 Y_0
!       .......................................
        USE VAPHFB_PAR
        implicit none
        integer Lam,jja,jla,jta,jjb,jlb,jtb
        real*8 delta_isos,fact,b2,radial,reduce
        double precision rnla
        REAL(DP) :: Wigner_3j
        integer jna,jnb,lja,ljb,tnlja,tnljb       
        integer    SPB_twoj_lj,SPB_l_lj

        !print *, HO%b_osc,FSPB%tnlj_fmax 
        if(.NOT. ALLOCATED(rME1B%T)) allocate(rME1B%T(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2))
        rME1B%T(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2) = 0.d0

           b2=HO%b_osc**2
        do Lam =0,2,2
           write(129,'(a40,i4)') 'The <tnlj1 ||Q_L|| tnlj2> with L=',Lam
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

          ! For both E0 and E2 transition: r^2 
          !radial=b2*rnla(jna,jla,Lam,jnb,jlb) !Rnl(jnla,jnlb)
          radial=b2*rnla(jna,jla,2,jnb,jlb) !Rnl(jnla,jnlb)

          if(Lam.eq.0 .and.jja.eq.jjb) rME1B%T(tnlja,tnljb,Lam) = radial*sqrt(jja+1.0) 
          if(Lam.eq.2) then
          fact = (iv(jla+jlb+Lam)+1.d0)/2.d0  
          rME1B%T(tnlja,tnljb,Lam)=radial*sqrt((2*Lam+1)/(4.*pi))&
      &                     *fact*iv((jjb+2*Lam-1)/2)              &
      &                     *sqrt((jja+1.)*(jjb+1.))          &
      &                     *Wigner_3j(jja,jjb,2*Lam,1,-1,0)
          endif

          if(abs(rME1B%T(tnlja,tnljb,Lam)).gt. CHOP) &
          write(129,'(3i6,f15.8)') tnlja,tnljb,Lam,rME1B%T(tnlja,tnljb,Lam)
        end do
        end do
        enddo ! Lam
        return

        end


       subroutine Calc_BE2()
       use VAPHFB_Par
       implicit none
       real*8 TD1B,fact,qpred,charge
       integer Ki,Kf,tnlja,tnljb,jta,jtb,jja,Ji,Jf,li,lf,Lam
       dimension TD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  !
       dimension qpred(0:JJmax,0:JJmax,0:2)
        integer jna,jnb,lja,ljb
        integer SPB_twoj_lj,SPB_l_lj

        Lam = 2
        print *, ' computing transitions between nuclear states'
        do Ji=GCM%Jmin,GCM%Jmax,GCM%Jdf !0,JJmax 
           ki = 0
        do Jf=GCM%Jmin,GCM%Jmax,GCM%Jdf !0,JJmax 
           kf = 0
           qpred(Ji,Jf,2) = 0.d0 
        do tnlja=1,FSPB%tnlj_fmax
           jta = FSPB%t(tnlja)  ! 0, 1
           jna = FSPB%n(tnlja) ! starting from 1 
           lja = FSPB%lj(tnlja)
           jja = SPB_twoj_lj(lja)    ! doubled
           fact = 1.0/sqrt(2*Lam+1.0)
        do tnljb=1,FSPB%tnlj_fmax

           jtb = FSPB%t(tnljb)

           if(jta.ne.jtb) cycle ! does not mix n with p

           if(Input%IntType.eq.1) then
              if(jta.eq.0) charge = 0.d0 
              if(jta.eq.1) charge = 1.d0
           else
              if(jta.eq.0) charge = Input%en
              if(jta.eq.1) charge = Input%ep
           endif

           qpred(Ji,Jf,2) = qpred(Ji,Jf,2) &
                                 + charge*fact *GCMDens%TD1BSum(tnlja,tnljb,Ji,Jf,1,1,2) &
                                 * rME1B%T(tnlja,tnljb,2)


        enddo 
        enddo 
        enddo 
        enddo 

        print *, '<J=2,1||Q2||J=2,1>=', qpred(2,2,2)
        print *, '<J=4,1||Q2||J=4,1>=', qpred(4,4,2)
        print *, '<J=6,1||Q2||J=6,1>=', qpred(6,6,2)
        print *, '<J=2,1||Q2||J=0,1>=', qpred(0,2,2) 
        print *, '<J=4,1||Q2||J=2,1>=', qpred(2,4,2) 
        print *, '<J=6,1||Q2||J=4,1>=', qpred(4,6,2) 

        return
        end
