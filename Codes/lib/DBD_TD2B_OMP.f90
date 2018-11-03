
        subroutine ZTD2B_Angles(ZRO,ZKAPA10,ZKAPA01,znorm,ZTBDJJ,NLEV)

        USE VAPHFB_PAR
        implicit none
        integer   NLEV 
        integer   k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 ztd,ZRO,ZKAPA10,ZKAPA01,ZTBDJJ,Znorm
        DIMENSION ZKAPA10(NLEV,NLEV)    ! density matrix
        DIMENSION ZKAPA01(NLEV,NLEV)    ! density matrix
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION ZTBDJJ    (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)

        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER a_ij,b_ij,a_kl,b_kl
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
        INTEGER b1,b2,bb1,bb2,ppnn
        INTEGER tt1,PP1,JJ1,tt2,PP2,JJ2

        INTEGER nlj1,nlj2,nlj3,nlj4
        REAL(DP) CG
        real*8 cg1,cg2

      do ppnn=1, ppnn_Max
!     ................................................ block bb1 for pp
           bb1= TPBO_ppnn%bb1(ppnn)
           b1 = TPBO_ppnn%PNonZero(bb1)
           tt1= TPB%btt(b1)
           PP1= TPB%bPP(b1)
           JJ1= TPB%bJJ(b1)
!      ................................................ check bb1=(tt,PP,JJ)
           if(bb1.ne.JJ1+PP1*(TPB%bJMax+1)+tt1*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'

!     ................................................ block bb2 for nn
           bb2= TPBO_ppnn%bb2(ppnn)
           b2 = TPBO_ppnn%NNonZero(bb2)
           tt2= TPB%btt(b2)
           PP2= TPB%bPP(b2)
           JJ2= TPB%bJJ(b2)

!      ................................................ check bb2
           if(bb2.ne.(JJ2+PP2*(TPB%bJMax+1)+tt2*2*(TPB%bJMax+1))) &
     &     stop ' Please check the basis'
!     .................................................
          if(JJ1 .ne. JJ2 .or.  PP1 .ne. PP2 .or. tt1.ne.2 .or. tt2.ne.0) &
     &    stop 'Errors in the basis for TD2B ..........'
!     .................................................
        JJ = JJ1
!-------------------------------------------------------
        do b_ij=0, TPB%aMax(bb1)-1
           a_ij= b_ij +1
           t_i = TPB%t1(bb1,a_ij)
           n_i = TPB%n1(bb1,a_ij)  ! n starts from 1 
           l_i = TPB%l1(bb1,a_ij)  ! l_i is doubled
           j_i = TPB%j1(bb1,a_ij)  ! j_i is doubled
           LJ_i = (l_i+j_i-1)/2

           t_j = TPB%t2(bb1,a_ij)
           n_j = TPB%n2(bb1,a_ij)  ! 
           l_j = TPB%l2(bb1,a_ij)  ! l_j is doubled
           j_j = TPB%j2(bb1,a_ij)  ! j_j is doubled
           LJ_j =(l_j+j_j-1)/2

           if(t_i+t_j .ne. 2) stop 'Error in the basis for TD2B'

        do b_kl=0, TPB%aMax(bb2)-1
           a_kl= b_kl  +1
           t_k = TPB%t1(bb2,a_kl)
           n_k = TPB%n1(bb2,a_kl)  ! 
           l_k = TPB%l1(bb2,a_kl)  ! l_i is doubled
           j_k = TPB%j1(bb2,a_kl)  ! j_i is doubled
           LJ_k = (l_k+j_k-1)/2

           t_l = TPB%t2(bb2,a_kl)
           n_l = TPB%n2(bb2,a_kl)  ! 
           l_l = TPB%l2(bb2,a_kl)  ! l_j is doubled
           j_l = TPB%j2(bb2,a_kl)  ! j_j is doubled
           LJ_l = (l_l+j_l-1)/2
           if(t_k+t_l .ne. 0) stop 'Error in the basis for TD2B'

!     ............... determine the nlj index
          nlj1 = SPB%nlj(n_i-1,lj_i)   ! n in the SPB start from 0
          nlj2 = SPB%nlj(n_j-1,lj_j)   ! while the n in ME2B_Base starts from 1
          nlj3 = SPB%nlj(n_k-1,lj_k)
          nlj4 = SPB%nlj(n_l-1,lj_l)

!      .......... the n in SPB%n starts from 0, while n in others starts from 1 !!!

       if(SPB%VType(lj_i,n_i-1,t_i) .ne. 1) cycle
       if(SPB%VType(lj_j,n_j-1,t_j) .ne. 1) cycle
       if(SPB%VType(lj_k,n_k-1,t_k) .ne. 1) cycle
       if(SPB%VType(lj_l,n_l-1,t_l) .ne. 1) cycle

!      .............. loops over m
       do m_i=-j_i,j_i,2
          k_i = tnljm%level(m_i,lj_i,n_i,t_i)
       do m_j=-j_j,j_j,2
          k_j = tnljm%level(m_j,lj_j,n_j,t_j)
          MM  = (m_i+m_j)/2

          cg1 = CG(j_i,m_i,j_j,m_j,2*JJ,2*MM)

       do m_k = -j_k,j_k,2
          m_l = MM*2 - m_k


!         do m_l=-j_l,j_l,2
!          if((m_k+m_l)/2 .ne. MM) cycle  ! A bug (not divided by 2) is found here
!         .......................................
          k_k = tnljm%level(m_k,lj_k,n_k,t_k)
          k_l = tnljm%level(m_l,lj_l,n_l,t_l)
          cg2=CG(j_k,m_k,j_l,m_l,2*JJ,2*MM)

!         .................................. cg2 = CG_Save((j_k+1)/2,(m_k+1)/2,(j_l+1)/2,(m_l+1)/2,JJ,MM) 

          ztd =  cg1*cg2*znorm/((2.d0*JJ+1.d0)*PNP%NFOM**2)       &
     &         *( ZRO(k_k,k_i)*ZRO(k_l,k_j)                   &
     &          -ZRO(k_l,k_i)*ZRO(k_k,k_j)                   &
     &          +ZKAPA01(k_i,k_j)*ZKAPA10(k_k,k_l))

          ZTBDJJ(a_ij,a_kl,ppnn) = ZTBDJJ(a_ij,a_kl,ppnn) + ztd

         enddo
         enddo
         enddo
!    ..................... end loops over m        
!          if(a_ij.eq.106 .and. a_ij .eq. a_kl .and.  ppnn.eq.3) then
         !   write(*,'(4i3,2i4,i3,f12.8)') nlj1,nlj2,nlj3,nlj4,a_ij,a_kl,JJ,dreal(ZTBDJJ(a_ij,a_kl,ppnn))
!            write(*,'(4i3,i3,f12.8)') nlj1,nlj2,nlj3,nlj4,JJ,dreal(ZTBDJJ(a_ij,a_kl,ppnn))
!          endif

        enddo
        enddo
        enddo ! ppnn
!    .................
        end


        subroutine ZTD2B_Integral(iangle,ZTBDJJ,J,zff,ZRho2BJJ)
        USE VAPHFB_PAR
        implicit none
        integer   iangle,J
        integer   bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        INTEGER b1,b2,bb1,bb2,ppnn
        INTEGER tt1,PP1,JJ1,tt2,PP2,JJ2
        complex*16 zff,ZRho2BJJ,ZTBDJJ
        DIMENSION ZTBDJJ    (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        DIMENSION ZRho2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)


      do ppnn=1, ppnn_Max
           bb1 = TPBO_ppnn%bb1(ppnn)
           b1  = TPBO_ppnn%PNonZero(bb1)
           tt1= TPB%btt(b1)
           PP1= TPB%bPP(b1)
           JJ1= TPB%bJJ(b1)
           if(bb1.ne.JJ1+PP1*(TPB%bJMax+1)+tt1*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'
!    ................................
           bb2 = TPBO_ppnn%bb2(ppnn)
           b2  = TPBO_ppnn%NNonZero(bb2)
           tt2= TPB%btt(b2)
           PP2= TPB%bPP(b2)
           JJ2= TPB%bJJ(b2)
           if(bb2.ne.JJ2+PP2*(TPB%bJMax+1)+tt2*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'
!    ................................
        do b_ij=0, TPB%aMax(bb1)-1
           a_ij= b_ij +1
        do b_kl=0, TPB%aMax(bb2)-1
           a_kl= b_kl  +1

!         attention: the "J" in weight(J) is
!         not the JJ in the two-body transition density
          ZRho2BJJ(a_ij,a_kl,ppnn) = ZRho2BJJ(a_ij,a_kl,ppnn)         &
     &                              + weightJ(J)*zff*ZTBDJJ(a_ij,a_kl,ppnn) 

         enddo
         enddo
         enddo
       end


        subroutine Write_ME2B_TD(iq1,iq2,ie,znorm1,znorm0,ZTD2BJJ,V_AB_CD_J)
        USE VAPHFB_PAR
        implicit none
        real*8     betac1,betac2,p001,p002
        integer    ie,iq1,iq2
        integer    bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        INTEGER    b1,b2,bb1,bb2,ppnn
        INTEGER    tt1,PP1,JJ1,tt2,PP2,JJ2
        real*8     rnorm,cNME(0:AMP%JJMax),cNME_tot
        real*8     d12,d34,factor
        integer    nlj1,nlj2,nlj3,nlj4        
        complex*16 zff,ZTD2BJJ,znorm0,znorm1,Rho2BN
        DIMENSION ZTD2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        real*8    V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax)
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
        complex*16 TBD
!     ...............................
      print *, 'JJmax=',HO%TwoJmax,AMP%JJMax
      open(300,file=File%TD2B,status='unknown')
 20   format(/)
 50   format('*** 2b DME from GCM calculation ***')
 60   format('# basis (block, aMax)')
195   format('bMax=',i6)
196   format(2i10)
197   format('# block',i6,':',' aMax=',i6)
198   format(2i10,f20.15)
199   format(7i5,2f12.7)
200   format(5i5,2f12.7)
250   format('*** Repeated 2B-DME:(',f12.7,') <-- (',2i6,')',6i3)
297   format('# pp_block',i6,':',' aMax=',i6,' # nn_block',&
     &       i6,':',' aMax=',i6) 

      write(300,*) dreal(znorm1),dreal(znorm0)  ! norm(q_i,q_i), norm(q_f,q_f)
      write(300,*)' a_ij   a_kl   bb  ppnn  TD2B'
      write(300,*)' -----------------------------------'

      write(ie,50)
      write(ie,60)
      write(ie,195) TPB%bMax-1
      if(TPB%bMax.gt.bMax) stop 'bMax is too small'
      do bm=0, TPB%bMax-1          ! The block starts from 0 to bMax_calc-1 
         write(ie,196) bm,TPB%cMax(bm) ! right 
      enddo

      write(ie,20)
!   ...................................
      cNME(0:AMP%JJMax) = zero
      do ppnn=1, ppnn_Max
           bb1 = TPBO_ppnn%bb1(ppnn)
           b1  = TPBO_ppnn%PNonZero(bb1)
           tt1= TPB%btt(b1)
           PP1= TPB%bPP(b1)
           JJ1= TPB%bJJ(b1)
           if(bb1.ne.JJ1+PP1*(TPB%bJMax+1)+tt1*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'
           bb2 = TPBO_ppnn%bb2(ppnn)
           b2  = TPBO_ppnn%NNonZero(bb2)
           tt2= TPB%btt(b2)
           PP2= TPB%bPP(b2)
           JJ2= TPB%bJJ(b2)
           if(bb2.ne.JJ2+PP2*(TPB%bJMax+1)+tt2*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'


         write(ie,297) b1,TPB%cMax(b1),b2,TPB%cMax(b2)

        do b_ij=0, TPB%aMax(bb1)-1
!-------------------------------------------------------
           a_ij= b_ij +1

        t_i = TPB%t1(bb1,a_ij)
        n_i = TPB%n1(bb1,a_ij)  ! 
        l_i = TPB%l1(bb1,a_ij)  ! l_i is doubled
        j_i = TPB%j1(bb1,a_ij)  ! j_i is doubled
        LJ_i = (l_i+j_i-1)/2

        t_j = TPB%t2(bb1,a_ij)
        n_j = TPB%n2(bb1,a_ij)  ! 
        l_j = TPB%l2(bb1,a_ij)  ! l_j is doubled
        j_j = TPB%j2(bb1,a_ij)  ! j_j is doubled
        LJ_j =(l_j+j_j-1)/2

        d12 = 0.d0
        if(t_i .eq. t_j .and. n_i .eq. n_j .and. LJ_i.eq.LJ_j) d12 = 1.d0

        do b_kl=0, TPB%aMax(bb2)-1
           a_kl= b_kl  +1
        t_k = TPB%t1(bb2,a_kl)
        n_k = TPB%n1(bb2,a_kl)  ! 
        l_k = TPB%l1(bb2,a_kl)  ! l_i is doubled
        j_k = TPB%j1(bb2,a_kl)  ! j_i is doubled
        LJ_k = (l_k+j_k-1)/2

        t_l = TPB%t2(bb2,a_kl)
        n_l = TPB%n2(bb2,a_kl)  ! 
        l_l = TPB%l2(bb2,a_kl)  ! l_j is doubled
        j_l = TPB%j2(bb2,a_kl)  ! j_j is doubled
        LJ_l = (l_l+j_l-1)/2



        d34 = 0.d0
        if(t_k .eq. t_l .and. n_k .eq. n_l .and. LJ_k.eq.LJ_l) d34 = 1.d0

        factor = dsqrt((1.d0+d12)*(1.d0+d34))
        rnorm  = dsqrt(dreal(znorm0*znorm1))
        Rho2BN = zzero
        if(abs(rnorm) .gt. 1.d-100)  Rho2BN = ZTD2BJJ(a_ij,a_kl,ppnn)/rnorm
        if(abs(rnorm) .lt. 1.d-100)  stop ' Norm is too small !'


!      .......... the n in SPB%n starts from 0, while n in others starts
!      from 1 !!!
!       A bug was here before: we should add the following four lines to
!       filter out the states outside the model space
!       if(max(nlj1,nlj2,nlj3,nlj4) .gt. HO%nljmax) print *,nlj1,nlj2,nlj3,nlj4,HO%nljmax
       if(SPB%VType(lj_i,n_i-1,t_i) .ne. 1) cycle
       if(SPB%VType(lj_j,n_j-1,t_j) .ne. 1) cycle
       if(SPB%VType(lj_k,n_k-1,t_k) .ne. 1) cycle
       if(SPB%VType(lj_l,n_l-1,t_l) .ne. 1) cycle

       nlj1 = SPB%nlj(n_i-1,lj_i)   ! n in the SPB start from 0
       nlj2 = SPB%nlj(n_j-1,lj_j)   ! while the n in ME2B_Base starts from 1
       nlj3 = SPB%nlj(n_k-1,lj_k)
       nlj4 = SPB%nlj(n_l-1,lj_l)
!      ............................................ calculate NME
       if(JJ1.le. min(HO%TwoJmax,AMP%JJMax)) then
!      V_AB_CD_J is the normalized ME2B
         cNME(JJ1) = cNME(JJ1) + 1.d0 &
     &                            *V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ1) &
     &                            *dreal(Rho2BN) *(2*JJ1+1.d0)/factor
       endif

!        if(a_ij.eq.38 .and. a_ij.eq.a_kl) &
!      & print *,V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ1),'Rho2B=',dreal(Rho2BN),factor

        if(abs(dreal(Rho2BN)).gt.1.d-7)&
     &  write(321,'(5i8,3f15.8)') nlj1,nlj2,nlj3,nlj4,JJ1, &
     &  sqrt(2*JJ1+1.d0)*dreal(Rho2BN),dreal(ZTD2BJJ(a_ij,a_kl,ppnn)),      &
     &  dreal(Rho2BN) 
 
!      ............................................ print out transition density
           if(ie.ne.300 .and. abs(Rho2BN).gt.CHOP) &
     &     write(ie,198) b_ij,b_kl,dreal(Rho2BN)

           if(abs(Rho2BN).gt.CHOP) &
     &     write(300,'(4i6,f15.10)')                      &
     &     a_ij,a_kl,bb2,ppnn,dreal(Rho2BN) !dreal(ZRho2BJJ(a_ij,a_kl,bb)/znorm)

        enddo  ! b_kl
        enddo  ! b_ij
        write(ie,20)
      enddo  ! ppnn
      
      close(300)
      close(ie)
      write(*,'(a,i2)') '....The 2B-TDME has been written to fort.',ie 

!     ......................... print out the J-decomposition of the NME
      cNME_tot = zero
      do JJ1=0,min(HO%TwoJmax,AMP%JJMax) !AMP%JJMax
         cNME_tot = cNME_tot + cNME(JJ1)
         write(*,'(a,i2,a,f12.8)') 'NME(J=',JJ1,')=',CNME(JJ1)
         write(31,'(a,i2,a,f12.8)') 'NME(J=',JJ1,')=',CNME(JJ1)
      enddo ! JJ1
         write(*,'(a,f12.8)') 'Total_NME=',cNME_tot
         write(31,'(a,f12.8)') 'Total_NME=',cNME_tot

          betac1 = Const1%beta2t_mesh(iq1)
          betac2 = Const%beta2t_mesh(iq2)
          p001  = Const1%P00_mesh(iq1)
          p002  = Const%P00_mesh(iq2)
         write(32,'(5f12.8)') betac1,p001,betac2,p002,cNME_tot

       end


!       ............................................
        subroutine DBD_NME_Calc(iq1,iq2,ZTD2BJJ,V_AB_CD_J)
        USE VAPHFB_PAR
        implicit none
        integer    ie,iq1,iq2
        integer    bm,tt,PP,JJ,bb,b_ij,a_ij,b_kl,a_kl
        INTEGER    b1,b2,bb1,bb2,ppnn
        INTEGER    tt1,PP1,JJ1,tt2,PP2,JJ2
        real*8     cNME(0:AMP%JJMax),cNME_tot
        real*8     d12,d34,factor,betac1,betac2,p001,p002
        integer    nlj1,nlj2,nlj3,nlj4
        complex*16 zff,ZTD2BJJ,ZRho2BN
        DIMENSION ZTD2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        real*8    V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax)
        INTEGER k_i,t_i,n_i,l_i,j_i,LJ_i,m_i
        INTEGER k_j,t_j,n_j,l_j,j_j,LJ_j,m_j
        INTEGER k_k,t_k,n_k,l_k,j_k,LJ_k,m_k
        INTEGER k_l,t_l,n_l,l_l,j_l,LJ_l,m_l
!     ...............................

      cNME(0:AMP%JJMax) = zero
      do ppnn=1, ppnn_Max
           bb1 = TPBO_ppnn%bb1(ppnn)
           b1  = TPBO_ppnn%PNonZero(bb1)
           tt1= TPB%btt(b1)
           PP1= TPB%bPP(b1)
           JJ1= TPB%bJJ(b1)
           if(bb1.ne.JJ1+PP1*(TPB%bJMax+1)+tt1*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'
           bb2 = TPBO_ppnn%bb2(ppnn)
           b2  = TPBO_ppnn%NNonZero(bb2)
           tt2= TPB%btt(b2)
           PP2= TPB%bPP(b2)
           JJ2= TPB%bJJ(b2)
           if(bb2.ne.JJ2+PP2*(TPB%bJMax+1)+tt2*2*(TPB%bJMax+1)) &
     &     stop ' Please check the basis'



        do b_ij=0, TPB%aMax(bb1)-1
!-------------------------------------------------------
           a_ij= b_ij +1

        t_i = TPB%t1(bb1,a_ij)
        n_i = TPB%n1(bb1,a_ij)  ! 
        l_i = TPB%l1(bb1,a_ij)  ! l_i is doubled
        j_i = TPB%j1(bb1,a_ij)  ! j_i is doubled
        LJ_i = (l_i+j_i-1)/2

        t_j = TPB%t2(bb1,a_ij)
        n_j = TPB%n2(bb1,a_ij)  ! 
        l_j = TPB%l2(bb1,a_ij)  ! l_j is doubled
        j_j = TPB%j2(bb1,a_ij)  ! j_j is doubled
        LJ_j =(l_j+j_j-1)/2

        d12 = 0.d0
        if(t_i .eq. t_j .and. n_i .eq. n_j .and. LJ_i.eq.LJ_j) d12 = 1.d0
        do b_kl=0, TPB%aMax(bb2)-1
           a_kl= b_kl  +1

        t_k = TPB%t1(bb2,a_kl)
        n_k = TPB%n1(bb2,a_kl)  ! 
        l_k = TPB%l1(bb2,a_kl)  ! l_i is doubled
        j_k = TPB%j1(bb2,a_kl)  ! j_i is doubled
        LJ_k = (l_k+j_k-1)/2

        t_l = TPB%t2(bb2,a_kl)
        n_l = TPB%n2(bb2,a_kl)  ! 
        l_l = TPB%l2(bb2,a_kl)  ! l_j is doubled
        j_l = TPB%j2(bb2,a_kl)  ! j_j is doubled
        LJ_l = (l_l+j_l-1)/2



       ZRho2BN = ZTD2BJJ(a_ij,a_kl,ppnn)


      d34 = 0.d0
      if(t_k .eq. t_l .and. n_k .eq. n_l .and. LJ_k.eq.LJ_l) d34 = 1.d0

       factor = dsqrt((1.d0+d12)*(1.d0+d34))
       nlj1 = SPB%nlj(n_i-1,lj_i)   ! n in the SPB start from 0
       nlj2 = SPB%nlj(n_j-1,lj_j)   ! while the n in ME2B_Base starts from 1
       nlj3 = SPB%nlj(n_k-1,lj_k)
       nlj4 = SPB%nlj(n_l-1,lj_l)
!      ............................................ calculate NME
         cNME(JJ1) = cNME(JJ1) +  V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ1) &
     &                            *dreal(ZRho2BN) *(2*JJ1+1.d0)/factor

        enddo  ! b_kl
        enddo  ! b_ij
      enddo  ! ppnn
!     ......................... print out the J-decomposition of the NME
      cNME_tot = zero
      do JJ1=0,AMP%JJMax
         cNME_tot = cNME_tot + cNME(JJ1)
         write(*,'(a,i2,a,f12.8)') 'NME(J=',JJ1,')=',CNME(JJ1)
         write(31,'(a,i2,a,f12.8)') 'NME(J=',JJ1,')=',CNME(JJ1)
      enddo ! JJ1
         write(*,'(a,f12.8)') 'Total_NME=',cNME_tot
         write(31,'(a,f12.8)') 'Total_NME=',cNME_tot

          betac1 = Const1%beta2t_mesh(iq1)
          betac2 = Const%beta2t_mesh(iq2)
          p001  = Const1%P00_mesh(iq1)
          p002  = Const%P00_mesh(iq2)
         write(32,'(5f12.8)') betac1,p001,betac2,p002,cNME_tot
198   format(2i10,2f12.8)

       end

        subroutine Read_ME2B(ie,ZRho2BJJ)
        USE VAPHFB_PAR
        implicit none
        integer     ppnn,ie,a12,a34,bb
        Complex*16  ZRho2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        real*8    cme

        ZRho2BJJ=zero 

        open(ie,file=File%TD2B,status='old')
        read(ie,*)
        read(ie,*)
        read(ie,*)
11      read(ie,*,end=21) a12,a34,bb,ppnn,cme
        ZRho2BJJ(a12,a34,ppnn) = cme
        goto 11
21      continue
        end
