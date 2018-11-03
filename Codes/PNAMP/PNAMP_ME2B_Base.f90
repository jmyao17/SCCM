      SUBROUTINE ME2B_Base_Full()
! ..............................................................................
!   New: added on March, 22, 2016, by jmyao
! ..............................................................................
!   This subroutine is to count the Two-Particle States
!   |(nlj)_1 (nlj)_2; J,PP,tt> in the block b=(J, PP,tt)
!   ...............
!   eMax=2*(n-1)+L
!   tt=0, 1, 2
!   PP=0, 1
!   n=1, 2, ..
!   L=0, 1, 2, Lmax
! ..............................................................................
      USE VAPHFB_Par 
      IMPLICIT NONE

      integer  bMax_tmp  
      INTEGER  JMax_tp(0:2,0:1)  
      INTEGER  a, bb, c, bm, J, PP, tt, total,me
      INTEGER  t1, n1, l1, j1, lj1, id1 
      INTEGER  t2, n2, l2, j2, lj2, id2
      INTEGER  J12m, J12p
      INTEGER  LJ_1, LJ_2
      INTEGER  Maxdim, n3,n4, n_1234
      INTEGER  idx(HO%nMax**4,1:4)
      INTEGER  icspb,tnlj_1,tnlj_2
!      INTEGER aMax(0:bMax)
!      INTEGER itpj 
!      itpj(tt,PP,J) = tt+PP*3+2*3*J
!     .............................. initialization
      TPB%aMax(0:bMax)            = 0
      JMax_tp(0:2,0:1)        = 0
      Maxdim                  = 0
      TPB%J12(0:bMax,1:aMaxMax)  = -1
      TPB%P12(0:bMax,1:aMaxMax) = -1
      TPB%n1(0:bMax,1:aMaxMax) = -1
      TPB%t1(0:bMax,1:aMaxMax) = -1
      TPB%twol1(0:bMax,1:aMaxMax) = -1
      TPB%j1(0:bMax,1:aMaxMax) = -1
      TPB%n2(0:bMax,1:aMaxMax) = -1
      TPB%t2(0:bMax,1:aMaxMax) = -1
      TPB%twol2(0:bMax,1:aMaxMax) = -1
      TPB%j2(0:bMax,1:aMaxMax) = -1
      TPB%iV(0:bMax,1:aMaxMax) = -1
      if(.NOT. allocated(TPB%block)) allocate(TPB%block(0:HO%twoJMax,0:1,0:2))
      TPB%block(0:HO%twoJMax,0:1,0:2) = -1  ! JJ, PP, TT
!     ............................. idx
       n_1234 = 0
       do n1=1,HO%nMax
       do n2=1,HO%nMax
       do n3=1,HO%nMax
       do n4=1,HO%nMax
          n_1234=n_1234 + 1
          idx(n_1234,1) = n1
          idx(n_1234,2) = n2
          idx(n_1234,3) = n3
          idx(n_1234,4) = n4
       enddo !
       enddo !
       enddo !
       enddo !
!     ............................. Full space for SPB
       if(.NOT. allocated(FSPB%tnlj)) allocate(FSPB%tnlj(0:HO%ljmax,0:HO%nmax,0:1))
       FSPB%tnlj(0:HO%ljmax,0:HO%nmax,0:1) = -1 
       icspb = 0 
       do t1=0,1
       do L1=0, HO%LMax*2, 2                ! L1 is doubled  
       do J1=abs(L1-1), (L1+1), 2     ! J1 is doubled
          LJ_1 = (L1+J1-1)/2
       do n1=1, (HO%eMax-L1/2)/2+1   ! eMax=2*(n1-1)+L1
          icspb = icspb + 1 

          FSPB%tnlj(LJ_1,n1,t1) = icspb 

       enddo
       enddo
       enddo
       enddo
       FSPB%tnlj_fmax = icspb
       if(.NOT. allocated(FSPB%t))     allocate(FSPB%t(1:FSPB%tnlj_fmax))
       if(.NOT. allocated(FSPB%n))     allocate(FSPB%n(1:FSPB%tnlj_fmax))
       if(.NOT. allocated(FSPB%lj))    allocate(FSPB%lj(1:FSPB%tnlj_fmax))
       icspb = 0
       do t1=0,1
       do L1=0, HO%LMax*2, 2                ! L1 is doubled  
       do J1=abs(L1-1), (L1+1), 2     ! J1 is doubled
          LJ_1 = (L1+J1-1)/2
       do n1=1, (HO%eMax-L1/2)/2+1   ! eMax=2*(n1-1)+L1
          icspb = icspb + 1

          FSPB%t(icspb) = t1
          FSPB%n(icspb) = n1
          FSPB%lj(icspb) = LJ_1

       enddo
       enddo
       enddo
       enddo



       if(.NOT. ALLOCATED(TPB%a)) ALLOCATE(TPB%a(0:bMax,0:icspb,0:icspb)) 
       TPB%a(0:bMax,0:icspb,0:icspb) = -1 

!      write(*,*) 'Construct TPB in Full Space ...'  
!     .............................. loop over ispsin (tt=t1+t2) 
      bMax_tmp = 0
      do tt=0, 2 ! nn, np, pp 
         t1 = int(tt/2)
         t2 = int((tt+1)/2)
      if(HO%lMax.gt.HO%eMax) stop 'ERROR: lMax is larger than eMax !'
      do L1=0, HO%LMax*2, 2                ! L1 is doubled  
      do J1=abs(L1-1), (L1+1), 2     ! J1 is doubled
         LJ_1 = (L1+J1-1)/2
        
         do L2=0, 2*HO%lMax, 2             ! L1 is doubled  
         do J2=abs(L2-1), (L2+1), 2     ! J1 is doubled
            LJ_2 = (L2+J2-1)/2
!---------------------------------------------
             J12m = abs(J1-J2)/2         ! J12m is not doubled
             J12p = abs(J1+J2)/2         ! J12p is not doubled
             PP   = mod((L1+L2)/2,2)     ! 0: positive; 1: negative
          do J=J12m, J12p            ! J is not doubled
             if(J+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1).gt.bMax_tmp) &
             bMax_tmp = J+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
!---------------------------------------------
            do n1=1, (HO%eMax-L1/2)/2+1   ! eMax=2*(n1-1)+L1
               id1 = n1 + (HO%NMax)*LJ_1
               tnlj_1 = FSPB%tnlj(LJ_1,n1,t1)
            do n2=1, (HO%eMax-L2/2)/2+1   ! eMax=2*(n2-1)+L2
               id2 = n2 + (HO%NMax)*LJ_2
               tnlj_2 = FSPB%tnlj(LJ_2,n2,t2)
!   ........................................ in the basis |ij>, only i<=j case 
            if(tt.ne.1.and.id1.gt.id2) cycle  ! only id1<=id2 is stored for nn and pp 
!         ....... two neutrons or two protons have the same quantum numbers (nlj), their
!                 angular momeum should be coupled to J=even in order to fulfill
!                 the anti-symmetry requirment for two identical particles.
             if(mod(J,2).eq.1.and.L1.eq.L2.and.J1.eq.J2.and.n1.eq.n2.and.tt.ne.1) cycle 
             if(J.gt.JMax_tp(tt,PP)) JMax_tp(tt,PP)= J
             bb = J+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
             if(bb.gt.bMax) stop 'bMax is too small' 
             TPB%aMax(bb) = TPB%aMax(bb) + 1
             a        = TPB%aMax(bb)
             if(a.gt.Maxdim) Maxdim = a            ! finding out the largest dimension of the blocks
             if(a.gt.aMaxMax) stop 'aMaxMax is too small'           
             if(bb.gt.bMax)    stop 'bMax    is too small'             
!    ------------------------------ b=0,...,; a=1,2.., 
              TPB%J12(bb,a)    = J                    
              TPB%P12(bb,a)    = PP       
              TPB%n1(bb,a)     = n1     ! n=1,2,  
              TPB%n2(bb,a)     = n2          
!              TPB%l1(bb,a)     = l1    ! doubled    
              TPB%twol1(bb,a)  = l1    ! doubled    
              TPB%j1(bb,a)     = j1    ! doubled   
              TPB%t1(bb,a)     = t1       
!              TPB%l2(bb,a)     = l2    ! doubled    
              TPB%twol2(bb,a)  = l2    ! doubled   
              TPB%j2(bb,a)     = j2    ! doubled   
              TPB%t2(bb,a)     = t2       

              TPB%a(bb,tnlj_1,tnlj_2) = a
              TPB%block(J,PP,tt)     = bb

            enddo
            enddo
          enddo ! J
          enddo ! J2
          enddo ! L2
        enddo ! J1
        enddo ! L1
       enddo ! tt
!     ............. print out
      bm = 0             ! counting the number of blocks
      total =0
      me    =0
      write(100,30)   
      write(100,40)   
      do tt=0,2
      do PP=0,1
      do J=0,JMax_tp(tt,PP) 
         bb = J+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
        if(TPB%aMax(bb).eq.0) cycle 
             TPB%btt(bm)   = tt
             TPB%bpp(bm)   = PP
             TPB%bJJ(bm)   = J
         write(100,10) bm, tt, PP, J,  TPB%aMax(bb)  ! bm=0,1,2,.. 
         TPB%cMax(bm) = TPB%aMax(bb)-1
         total = total + TPB%aMax(bb)
         me    = me + TPB%aMax(bb)**2
        
!    ......................................
!         if(tt.eq.0.and.PP.eq.0.and.J.eq.0) then
         if(bm.eq.27) then
          write(101,*) 'aMax=',TPB%cMax(bm)
          write(101,*) 'a    n1   l1   j1  n2  l2  j2   JJ  PP  t1   t2'
          do a=1, TPB%aMax(bb)  ! loop over the two-particle states in block b 
             write(101,5) a-1,TPB%n1(bb,a),TPB%twol1(bb,a)/2,TPB%j1(bb,a),&
      &                   TPB%n2(bb,a),TPB%twol2(bb,a)/2,TPB%j2(bb,a),&
      &         TPB%J12(bb,a),TPB%P12(bb,a),TPB%t1(bb,a),TPB%t2(bb,a)   
          enddo ! a
         endif
!    ......................................
         bm =bm + 1
       enddo ! J
       enddo ! PP 
       enddo ! tt

       TPB%bMax = bm
       if(TPB%bMax.gt.bMax) stop 'bMax for ME2B is too small !'
       write(100,40)   
       write(100,20) bm-1, Maxdim, total, me  
       write(100,*) 'bMax (including temp. max.value):',max(bm,bMax_tmp)
       write(100,40)   
  1    format('(  n1 l1  j1,  n2 l2  j2) ->  J  Pi |t1  t2')
  5    format('a=',i3,'(',3i3,'/2,',3i3,'/2',') ->',2i3,'  |',2i3)
 10    format(i6,':',3i6,'  |',i6)
 20    format('b= 0,...,',i3,'|  Maxdim=',i3,'| total dim:',i8,&
        &     '| No. of me:',i12)
 30    format('     b:     tt    Pi    J|   dim')
 40    format('------------------------------------------------------')
!      .......
       Return
       END SUBROUTINE 
!   ....................................................
       Integer function TPBC_bCalc(tt, PP, J, JMax_tp)
       Implicit none
       INTEGER tt, PP, J, JMax_tp
       TPBC_bCalc=((tt*2 + PP)*(JMax_tp+1) + J) 
       return
       End
