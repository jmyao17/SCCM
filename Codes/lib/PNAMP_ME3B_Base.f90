!    ************************************************
!    Construct a basis for three-body matrix elements
!    in a given model space
!    ATTENTION:
!    n=0,1,2,3,...
!    ************************************************
           subroutine Rho3B_1B_Model_Space(lpr)
!    ...................................................
!          One-body basis
!    ...................................................
           USE VAPHFB_PAR
           implicit none

           logical lpr
           integer a,nlj,nlj_vmax,tnlj_vmax,lj,tnlj,t,n,l,twoj
           integer NOrbit(0:1)
           integer a_1,t_1,L_1,J_1,N_1,LJ_1,tnlj_1

        if(lpr) then
          write(*,*) '*******************************************' 
          write(*,*) ' The basis for the Rho3B in a Valence Space' 
          write(*,*) '*******************************************' 
           write(*,*) '--------------------------------'
           write(*,*) ' Definition of Valence Sapce    '
           write(*,*) '--------------------------------'
           write(*,*) '     t,   n,   l,  2j,  2j+1'
         endif
           NOrbit(0) = 0
           NOrbit(1) = 0
!        ........................ read model (valence) space
           open(31, file="../val/me3b_"//trim(Input%vs4me3b)//".val",status="old")
           read(31,*)
           read(31,*) tnlj_vmax
           print *, tnlj_vmax,HO%tnljmax
           Rho3B%VSPB%tnlj_vmax = tnlj_vmax
           if(Rho3B%VSPB%tnlj_vmax .le.0 .or. tnlj_vmax .gt. HO%tnljmax) &
     &     stop 'dimension is not properly defined !'
           if(.NOT. ALLOCATED(Rho3B%VSPB%t))  ALLOCATE(Rho3B%VSPB%t(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%n))  ALLOCATE(Rho3B%VSPB%n(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%l))  ALLOCATE(Rho3B%VSPB%l(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%twoj))  ALLOCATE(Rho3B%VSPB%twoj(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%lj))    ALLOCATE(Rho3B%VSPB%lj(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%VType)) ALLOCATE(Rho3B%VSPB%VType(1:tnlj_vmax))
           if(.NOT. ALLOCATED(Rho3B%VSPB%tnlj))  ALLOCATE(Rho3B%VSPB%tnlj(0:HO%ljmax,0:HO%nmax,0:1))
!        ...................... initialization
          do a=1,tnlj_vmax
              Rho3B%VSPB%t(a)     = -1
              Rho3B%VSPB%n(a)     = -1
              Rho3B%VSPB%l(a)     = -1
              Rho3B%VSPB%twoj(a)  = -1
              Rho3B%VSPB%lj(a)    = -1
              Rho3B%VSPB%VType(a) = -1
          enddo
           read(31,*)
  1        read(31,*,end=2) tnlj, t, n, l, twoj

           if(tnlj.gt.tnlj_vmax) stop 'me3j.val is defined improperly !'
              Rho3B%VSPB%t(tnlj)     = t    ! t=0(n),1(p)
              Rho3B%VSPB%n(tnlj)     = n    ! n=0,1,2,..
              Rho3B%VSPB%l(tnlj)     = l    ! l is NOT doubled
              Rho3B%VSPB%twoj(tnlj)  = twoj
              lj = (l*2 + twoj -1)/2
              Rho3B%VSPB%tnlj(lj,n,t)= tnlj
              Rho3B%VSPB%lj(tnlj)    = lj
              Rho3B%VSPB%VType(tnlj) = 1
              if(t.eq.0) NOrbit(0) = NOrbit(0) + (Rho3B%VSPB%twoj(tnlj)+1)
              if(t.eq.1) NOrbit(1) = NOrbit(1) + (Rho3B%VSPB%twoj(tnlj)+1)
              write(*,'(5i6)') t, Rho3B%VSPB%n(tnlj),Rho3B%VSPB%l(tnlj), &
     &                        Rho3B%VSPB%twoj(tnlj), Rho3B%VSPB%twoj(tnlj)+1
              goto 1
  2        continue
         if(lpr) write(*,*) '--------------------------------'

         if(NOrbit(0).ne.NOrbit(1))  &
     &   stop 'Numbers of levels for Neutron and protons are not equal!'

         Rho3B%vnlj%nlj_vmax         =  tnlj_vmax/2
         nlj_vmax                    =  tnlj_vmax/2

         if(.NOT. ALLOCATED(Rho3B%vnlj%n))  ALLOCATE(Rho3B%vnlj%n(1:nlj_vmax))
         if(.NOT. ALLOCATED(Rho3B%vnlj%l))  ALLOCATE(Rho3B%vnlj%l(1:nlj_vmax))
         if(.NOT. ALLOCATED(Rho3B%vnlj%lj))  ALLOCATE(Rho3B%vnlj%lj(1:nlj_vmax))
         if(.NOT. ALLOCATED(Rho3B%vnlj%twoj))  ALLOCATE(Rho3B%vnlj%twoj(1:nlj_vmax))
         do tnlj=1,tnlj_vmax/2
              Rho3B%Vnlj%n(tnlj)      = Rho3B%VSPB%n(tnlj)
              Rho3B%Vnlj%l(tnlj)      = Rho3B%VSPB%l(tnlj)  ! L is NOT doubled
              Rho3B%Vnlj%lj(tnlj)     = Rho3B%VSPB%lj(tnlj)
              Rho3B%Vnlj%twoj(tnlj)   = Rho3B%VSPB%twoj(tnlj)
!              write(*,'(4i4)') tnlj,Rho3B%Vnlj%n(tnlj),Rho3B%Vnlj%l(tnlj),Rho3B%Vnlj%twoj(tnlj)
         enddo

       end  

!    ...................................................
           subroutine Rho3B_2B_Model_Space(lpr)
!    ...................................................
!          One-body basis
!    ...................................................
           USE VAPHFB_PAR
           implicit none

           logical lpr
           INTEGER  JMax_tp(0:2,0:1)
           INTEGER  a, bb, c, bm, JJ, PP, tt, total,me
           INTEGER  n_1, l_1, j_1, id1
           INTEGER  n_2, l_2, j_2, id2
           INTEGER  J12m, J12p
           INTEGER  LJ_1, LJ_2, LJ_3
           INTEGER  tnlj_1,tnlj_2
           INTEGER  Maxdim,n4,n_1234

           INTEGER t_1,t_2,a_1,a_2,bc
           INTEGER tnlj,sp_tnljm,spidx
           INTEGER t,n,l,lj,twoj,twom,tnljm_dMax
!    ...................................................
        Rho3B%vtnljm%tnljm_vmax = 0
        tnljm_dMax = 0 
        do tnlj=1, Rho3B%VSPB%tnlj_vmax  ! 8
            twoj = Rho3B%VSPB%twoj(tnlj)    ! doubled   
            do twom = -twoj,twoj,2  ! doubled  
               tnljm_dMax = tnljm_dMax +1
            enddo ! m_j3
       enddo ! a
       Rho3B%vtnljm%tnljm_vmax = tnljm_dmax


       if(lpr) then
          write(*,*) ' Number of S.P. States in Valence Space ' 
          write(*,'(a,i4)') ' Labeled with tnlj :',Rho3B%VSPB%tnlj_vmax
          write(*,'(a,i4)') ' Labeled with tnljm:',tnljm_dMax
          write(*,*) '---------------------------------------'
          write(*,*) ' idx,  t,  n,  l, 2j,  2m+1'
          write(*,*) '---------------------------------------'
       endif
      if(.NOT. ALLOCATED(Rho3B%vtnljm%t))  ALLOCATE(Rho3B%vtnljm%t(1:tnljm_dmax))
      if(.NOT. ALLOCATED(Rho3B%vtnljm%n))  ALLOCATE(Rho3B%vtnljm%n(1:tnljm_dmax))
      if(.NOT. ALLOCATED(Rho3B%vtnljm%l))  ALLOCATE(Rho3B%vtnljm%l(1:tnljm_dmax))
      if(.NOT. ALLOCATED(Rho3B%vtnljm%twoj))  ALLOCATE(Rho3B%vtnljm%twoj(1:tnljm_dmax))
      if(.NOT. ALLOCATED(Rho3B%vtnljm%twom))  ALLOCATE(Rho3B%vtnljm%twom(1:tnljm_dmax))
      if(.NOT. ALLOCATED(Rho3B%vtnljm%lj))  ALLOCATE(Rho3B%vtnljm%lj(1:tnljm_dmax))
       if(.NOT. ALLOCATED(Rho3B%vtnljm%tnljm))  &
     & ALLOCATE(Rho3B%vtnljm%tnljm(-HO%twojmax:HO%twojmax,0:HO%ljmax,0:HO%nmax,0:1))


        if(.NOT. ALLOCATED(Rho3B%VTPB%block)) &
     &  ALLOCATE(Rho3B%VTPB%block(0:HO%twojmax,0:1,0:2))
        if(.NOT. ALLOCATED(Rho3B%VTPB%a)) &
     &  ALLOCATE(Rho3B%VTPB%a(0:bMax,1:Rho3B%VSPB%tnlj_vmax,1:Rho3B%VSPB%tnlj_vmax))
        tnljm_dMax = 0 
        do tnlj=1, Rho3B%VSPB%tnlj_vmax  ! 8
            twoj = Rho3B%VSPB%twoj(tnlj)    ! doubled   
            do twom = -twoj,twoj,2  ! doubled  
               tnljm_dMax = tnljm_dMax +1
               t  = Rho3B%VSPB%t(tnlj)
               l  = Rho3B%VSPB%l(tnlj)
               n  = Rho3B%VSPB%n(tnlj)
               lj = Rho3B%VSPB%lj(tnlj)
               Rho3B%vtnljm%tnljm(twom,lj,n,t) = tnljm_dMax

               Rho3B%vtnljm%t(tnljm_dMax)  = t 
               Rho3B%vtnljm%n(tnljm_dMax)  = n 
               Rho3B%vtnljm%l(tnljm_dMax)  = l          ! L is NOT doubled 
               Rho3B%vtnljm%twoj(tnljm_dMax)  = twoj
               Rho3B%vtnljm%lj(tnljm_dMax)    = lj
               Rho3B%vtnljm%twom(tnljm_dMax)  = twom
               if(lpr)    write(*,'(6i4)') tnljm_dmax,t,n,l,twoj,twom

            enddo ! m_j3
       enddo ! a
!    ----------------------------------
!    Two-Particle Basis |(a_1 a_2)>
!    ----------------------------------
!     .............................. initialization
      Rho3B%VTPB%aMax(0:bMax) = 0
      JMax_tp(0:2,0:1)        = 0
      Maxdim                  = 0
      do tt=0, 2 ! nn, np, pp 
         t_1 = int(tt/2)
         t_2 = int((tt+1)/2)
         do a_1=1, Rho3B%Vnlj%nlj_vmax
            L_1  = Rho3B%Vnlj%l(a_1)*2   ! twol
            J_1  = Rho3B%Vnlj%twoj(a_1)  ! twoj
            N_1  = Rho3B%Vnlj%n(a_1)     ! starts from 0
            LJ_1 = Rho3B%Vnlj%lj(a_1)
            if(LJ_1 .ne. (L_1+J_1-1)/2) stop 'Error in ME3B_Base !'
            tnlj_1 = Rho3B%VSPB%tnlj(LJ_1,N_1,t_1)

            id1 = N_1 + (HO%nmax)*LJ_1
         do a_2=1,  Rho3B%Vnlj%nlj_vmax
            L_2  = Rho3B%Vnlj%l(a_2)*2   ! twol
            J_2  = Rho3B%Vnlj%twoj(a_2)  ! twoj
            N_2  = Rho3B%Vnlj%n(a_2)
            LJ_2 = Rho3B%Vnlj%lj(a_2)
            LJ_2 = (L_2+J_2-1)/2
!            tnlj_2 = LJ_2 + (N_2-1)*LJMax + t_2*NMax*LJMax
            tnlj_2 = Rho3B%VSPB%tnlj(LJ_2,N_2,t_2)
            id2 = N_2 + (HO%nmax)*LJ_2
! .... test
!           if(tt.ne.1.and.id1.gt.id2) cycle
! .... test
             J12m = abs(J_1-J_2)/2         ! J12m is not doubled
             J12p = abs(J_1+J_2)/2         ! J12p is not doubled
             PP   = mod((L_1+L_2)/2,2)     ! 0: positive; 1: negative
           do JJ=J12m, J12p            ! J is not doubled
              if(JJ.gt.JMax_tp(tt,PP)) JMax_tp(tt,PP)= JJ
  
             bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
 
            if(bb.gt.bMax) then
              write(*,*) 'Error in SM_Base: bMax should be larger than:',bb
              stop
             endif
!      .......... counter
             Rho3B%VTPB%aMax(bb) = Rho3B%VTPB%aMax(bb) + 1
!             write(*,'(6i5)') tnlj_1,tnlj_2,JJ,PP,bb,Rho3B%VTPB%aMax(bb)
             a        = Rho3B%VTPB%aMax(bb)
             if(a.gt.Maxdim) Maxdim = a            ! finding out the largest dimension of the blocks
             if(a.gt.aMaxMax)  stop 'aMaxMax is too small'
             if(bb.gt.bMax)    stop 'bMax    is too small'
!    ------------------------------ b=0,...,; a=1,2.., 
             Rho3B%VTPB%J12(bb,a)    = JJ     ! Not doubled
             Rho3B%VTPB%P12(bb,a)    = PP
             Rho3B%VTPB%n1(bb,a)     = N_1
             Rho3B%VTPB%n2(bb,a)     = N_2
             Rho3B%VTPB%twol1(bb,a)  = L_1    ! doubled    
             Rho3B%VTPB%j1(bb,a)     = J_1    ! doubled   
             Rho3B%VTPB%t1(bb,a)     = t_1
             Rho3B%VTPB%twol2(bb,a)  = L_2    ! doubled   
             Rho3B%VTPB%j2(bb,a)     = J_2    ! doubled   
             Rho3B%VTPB%t2(bb,a)     = t_2
             Rho3B%VTPB%a(bb,tnlj_1,tnlj_2) = a
             Rho3B%VTPB%block(JJ,PP,tt)     = bb
            if(JJ.gt. HO%twojmax) stop 'twojmax is too small'
            if(lpr) write(401,'(4i5)') JJ,PP,tt,Rho3B%VTPB%block(JJ,PP,tt)
          enddo
          enddo
      enddo
      enddo

!-----------------------------
! counting the # of the blocks
!-----------------------------
      bc    =0             ! counting the number of blocks
      total =0
      me    =0
      if(lpr) then
        write(200,30)
        write(200,40)
        write(*,*) '---------------------------------------'
        write(*,*) ' Definition of TPB in Valence Sapce    '
        write(*,*) '---------------------------------------'
        write(*,'(a)') '   tt   PP   JJ   bb   aMax(bb)'
        !write(*,*) 'bJMax=',TPB%bJMax
      endif
      do tt=0,2
      do PP=0,1
      do JJ=0,JMax_tp(tt,PP)
         bb = JJ+PP*(TPB%bJMax+1)+tt*2*(TPB%bJMax+1)
         if(Rho3B%VTPB%aMax(bb).eq.0) cycle
             Rho3B%VTPB%btt(bc)   = tt
             Rho3B%VTPB%bpp(bc)   = PP
             Rho3B%VTPB%bJJ(bc)   = JJ
         if(lpr) write(200,10) bc, tt, PP, JJ,  Rho3B%VTPB%aMax(bb)  ! bc=0,1,2,.. 
         if(lpr) write(*,'(5i5)') tt,PP,JJ,bb,Rho3B%VTPB%aMax(bb) !,tnlj_1,tnlj_2,pTPB%a(bb,tnlj_1,tnlj_2)
         Rho3B%VTPB%cMax(bc) = Rho3B%VTPB%aMax(bb)-1
         total = total + Rho3B%VTPB%aMax(bb)
         me    = me + Rho3B%VTPB%aMax(bb)**2

!    ......................................
!         if(tt.eq.0.and.PP.eq.0.and.JJ.eq.0) then
         if(lpr .and. bb.eq.2) then
           write(201,*) 'aMax=',Rho3B%VTPB%cMax(bc)
           write(201,*) '   a  n1 l1 j1  n2  l2  j2   JJ  PP  t1   t2'
          do a=1, Rho3B%VTPB%aMax(bb)  ! loop over the two-particle states in block b 
             write(201,5) a-1,Rho3B%VTPB%n1(bb,a),Rho3B%VTPB%twol1(bb,a)/2,Rho3B%VTPB%j1(bb,a),&
      &                   Rho3B%VTPB%n2(bb,a),Rho3B%VTPB%twol2(bb,a)/2,Rho3B%VTPB%j2(bb,a),&
      &         Rho3B%VTPB%J12(bb,a),Rho3B%VTPB%P12(bb,a),Rho3B%VTPB%t1(bb,a),Rho3B%VTPB%t2(bb,a)
          enddo ! a
         endif
!    ......................................
         bc =bc + 1
       enddo ! J
       enddo ! PP 
       enddo ! tt
       Rho3B%VTPB%bMax = bc
!    ......................................
       if(lpr) then
         write(200,40)
         write(200,20) bc-1, Maxdim, total, me
         write(200,40)
       endif
  1    format('(  n1 l1  j1,  n2 l2  j2) ->  J  Pi |t1  t2')
  5    format('a=',i3,'(',3i3,'/2,',3i3,'/2',') ->',2i3,'  |',2i3)
 10    format(i6,':',3i6,'  |',i6)
 20    format('b= 0,...,',i3,'|  Maxdim=',i3,'| total dim:',i8,&
        &     '| No. of me:',i12)
 30    format('     b:     tt    Pi    J|   dim')
 40    format('------------------------------------------------------')



       return
      END




      SUBROUTINE ME3B_Base_Full(lpr)
      USE VAPHFB_Par
      IMPLICIT NONE
      logical  lpr
      INTEGER  bmax_calc
      INTEGER  a, bb, c, bm, J, PP, tt, total,me
      INTEGER  t1, n1, l1, j1, lj1, id1
      INTEGER  t2, n2, l2, j2, lj2, id2
      INTEGER  t3, n3, l3, j3, lj3, id3
      INTEGER  J12m, J12p
      INTEGER  LJ_1, LJ_2, LJ_3
      INTEGER  m_j3
      INTEGER  a12,J12,P12,t12

      INTEGER J123,P123,t123,J12m3,J12p3
      INTEGER d123,a3,spidx,tnljm_dMax
      INTEGER nljm1,nljm2,nljm3
!     .............................. 
!        spidx = 0
!        do nljm1=1, p1Bm%tnljm_dMax  ! 
!        do nljm2=1, p1Bm%tnljm_dMax  ! 
!        do nljm3=1, p1Bm%tnljm_dMax  ! 
!           spidx = spidx +1
!           p3B_idx(nljm1,nljm2,nljm3) = spidx
!        enddo
!        enddo
!        enddo
!       write(*,'(a,i6,a,i6,a)') 'Dimension of (1,2,3) in m-scheme:', &
!      &     spidx,'(Max:',spMaxMax3,')'
!       if(spidx.gt.spMaxMax3)  &
!     & stop 'Error in SM_Base: spMaxMax3 should be larger'

!      ......................
!      Three-particle states 
!      .....................

          if(.NOT. allocated(Rho3B%idx)) allocate(Rho3B%idx(0:bMax,1:aMaxMax,1:HO%tnljmax,0:HO%twojmax*3))
          Rho3B%idx(0:bMax,1:aMaxMax,1:HO%tnljmax,0:HO%twojmax*3) = -1
          d123 = 0
          do bm=0, Rho3B%VTPB%bMax-1
             t12= Rho3B%VTPB%btt(bm)
             P12= Rho3B%VTPB%bPP(bm)
             J12= Rho3B%VTPB%bJJ(bm)  ! doubled
             bb = J12+P12*(TPB%bJMax+1)+t12*2*(TPB%bJMax+1)
             do a12=1,Rho3B%VTPB%aMax(bb),1
                do a3=1, Rho3B%VSPB%tnlj_vmax
                   J12m3 = abs(2*J12-Rho3B%VSPB%twoj(a3)) ! doubled
                   J12p3 = abs(2*J12+Rho3B%VSPB%twoj(a3)) ! doubled 
                do J123=J12m3, J12p3, 2            ! doubled 

                   d123 = d123 + 1
                enddo ! J123   
!-------------------------------------------------------
             enddo ! a12
       enddo ! tt
       enddo ! tt
        Rho3B%idx123_vmax = d123

      if(lpr) then
         write(*,*) '-------------------------------------------------'
         write(*,'(a,i6,a,i6,a)') ' Dimension of Rho3B (m-scheme):',d123
         write(*,*) '-------------------------------------------------'
         write(*,*) ' The basis for the Rho3B written into fort.400.'
         write(400,*) '             bm,  bb,  a12,  a3,  J123'
      endif
      if(.NOT. ALLOCATED(Rho3B%t12))  ALLOCATE(Rho3B%t12(1:d123))
      if(.NOT. ALLOCATED(Rho3B%J12))  ALLOCATE(Rho3B%J12(1:d123))
      if(.NOT. ALLOCATED(Rho3B%P12))  ALLOCATE(Rho3B%P12(1:d123))
      if(.NOT. ALLOCATED(Rho3B%a12))  ALLOCATE(Rho3B%a12(1:d123))
      if(.NOT. ALLOCATED(Rho3B%t3))  ALLOCATE(Rho3B%t3(1:d123))
      if(.NOT. ALLOCATED(Rho3B%n3))  ALLOCATE(Rho3B%n3(1:d123))
      if(.NOT. ALLOCATED(Rho3B%lj3))  ALLOCATE(Rho3B%lj3(1:d123))
      if(.NOT. ALLOCATED(Rho3B%J123))  ALLOCATE(Rho3B%J123(1:d123))

          d123 = 0
          do bm=0, Rho3B%VTPB%bMax-1
             t12= Rho3B%VTPB%btt(bm)
             P12= Rho3B%VTPB%bPP(bm)
             J12= Rho3B%VTPB%bJJ(bm)  ! doubled
             bb = J12+P12*(TPB%bJMax+1)+t12*2*(TPB%bJMax+1)
!                write(*,*) 'aMax(',bm,bb,t12,bJJ(bm),')=',aMax(bb)
             do a12=1,Rho3B%VTPB%aMax(bb),1

                do a3=1, Rho3B%VSPB%tnlj_vmax
                   n3 = Rho3B%VSPB%n(a3)    ! L is NOT double
                   l3 = Rho3B%VSPB%l(a3)    ! L is NOT double
                   j3 = Rho3B%VSPB%twoj(a3) ! double
                   LJ_3 = (l3*2+j3-1)/2 ! lj

                   J12m3 = abs(2*J12-Rho3B%VSPB%twoj(a3)) ! doubled
                   J12p3 = abs(2*J12+Rho3B%VSPB%twoj(a3)) ! doubled 
!                   P123  = P12*(-1)**(p1B%l(a3)/2)
                   t123  = t12+Rho3B%VSPB%t(a3)      !   
                do J123=J12m3, J12p3, 2            ! doubled 

                   d123 = d123 + 1
!                   if(d123.gt.dMaxMax) &
!     &             stop 'Error: dMaxMax should be larger!'
                   Rho3B%t12(d123) = t12
                   Rho3B%J12(d123) = J12             ! is not doubled 
                   Rho3B%P12(d123) = P12
                   Rho3B%a12(d123) = a12
                   Rho3B%t3(d123)  = Rho3B%VSPB%t(a3)
                   Rho3B%n3(d123)  = Rho3B%VSPB%n(a3)
                   Rho3B%lj3(d123) = LJ_3
                   Rho3B%J123(d123)= J123            ! doubled
                   Rho3B%idx(bb,a12,a3,J123)= d123
                   if(lpr) write(400,'(a,5i5,a,i8)') 'Rho3B%idx(',bm,bb,a12,a3,J123,')=',d123
!           if(d123.eq.2363) then
!             write(*,*) d123, 'J123=',J123,' t12=',t12,' J12=',J12
!            write(*,*) '(123):', t3,n3,l3,j3
!           endif
                enddo ! J123   
!-------------------------------------------------------
             enddo ! a12
       enddo ! tt
       enddo ! tt
 10   format(a,3i4,',',3i4,',',3i4)
       !stop 'ME3B_Base!!'
       return
       END  
