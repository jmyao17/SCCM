!          ................................................................
!          Read the ME2B files (in J-scheme) from Nathan and transform into
!          <(12)J;TMT | V | (34)J;TMT>
!          .................................................................
           subroutine EOM_ME2B_J2M()
           USE VAPHFB_PAR
           implicit none

        integer jja,jma,jta,jna,jlja,jacoup
        integer jjb,jmb,jtb,jnb,jljb,jbcoup
        integer jjc,jmc,jtc,jnc,jljc,jccoup
        integer jjd,jmd,jtd,jnd,jljd,jdcoup
        integer MJ,MMJ,jmax,nlj,iabcd,icheck
        real*8  AN_INV,AN_AB,AN_CD,cb1,cb2,cb3,cb4,suma,sumando
        real*8  ajtot,attot,ame_coup
        real*8  phasJT,phasab,phascd,delta_ab,delta_cd


           REAL(DP) CG
           integer tnlj,tz1
           integer t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ
           integer twoj1,twoj2,iphase_12   
           integer l1,l2,l3,l4
           integer j1,id1,id2,id3,id4
           integer twoj3,twoj4,iphase_34    
           real*8  fac_cg,me,cg_coeff(-1:1,-1:1,0:1),cb,occ
           real*8, dimension(:,:,:,:,:,:), allocatable :: V_AB_CD_J
           integer ia,ib,ic,id,iT,iJ,iphase_cd
           integer iphaseJ,iphaseJT
           integer nlj1,nlj2,nlj3,nlj4,MMT,MT,TT
           real*8  delta_J, delta_JT,d_ab,d_cd,d_abcd
           real*8  d_12,d_34,d_1234
           real*8  AN_12,AN_34
           integer it1234(0:1,0:1,0:1,0:1)
           integer itt1(0:5),itt2(0:5),itt3(0:5),itt4(0:5)
           integer ii,it1,it2,it3,it4
           integer iphase_CS
           integer iprint
           real*8 diff
           character(500) find_file

           integer, dimension(:), allocatable :: sps_t,sps_n,sps_lj

          INT_DIR = find_file("GCM_ME_FILES",File%IMSRG_Hme2b)
          File%f2bJ=File%IMSRG_Hme2b

!     .............. initialization
       if(.NOT. ALLOCATED(V_AB_CD_J))  &
     & ALLOCATE(V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5))

!      ............. proton(-1); neutron (+1)
          do TT= 0,HO%Tmax  ! 1  TT is NOT doubled
          do t2=-1,1,2        ! t2 is doubled
          do t1=-1,1,2        ! t1 is doubled
             cg_coeff(t1,t2,TT) = zero
             MMT = t1+t2          ! MMT is doubled
             cb = CG(1,t1,1,t2,2*TT,MMT)
             cg_coeff(t1,t2,TT) = cb
          enddo
          enddo
          enddo
         write(*,*) 'nlj   lj   twoj'

! .....................................
! 0: nnnn
! 1: npnp
! 2: nppn
! 3: pnnp
! 4: pnpn
! 5: pppp
! .....................................
       ii = 0
      do it1=0,1
      do it2=0,1
      do it3=0,1
      do it4=0,1
         if(it1+it2 .ne. it3+it4) cycle
         it1234(it1,it2,it3,it4) = ii
         itt1(ii) = it1
         itt2(ii) = it2
         itt3(ii) = it3
         itt4(ii) = it4
         print*,it1,it2,it3,it4,ii
         ii = ii + 1
      enddo
      enddo
      enddo
      enddo

!      only two-body
!       if(Input%IntType.eq.1) then
!          if(Input%InME3B.ne.1) 
!          if(Input%InME3B.eq.1) & 
!        &  File%f2bJ='../Int/IMSRG_'//Input%cFlow//'_chi2b00000_srg0625'//&
!        &  Input%ctpp//Input%cValID//'_'//Input%chwHO//'_J.dat'
!       else
!          f2bJ_save='../Int/IMSRG_'//Input%cFlow//'_'&
!     &    //Input%cIntID//'_J.dat'
!       endif
!        open(300,file=f2bJ_save,status='unknown')

        write(*,*) '    Read ME2B from ',trim(INT_DIR)//File%f2bJ
        open(26,file=trim(INT_DIR)//File%f2bJ,status='old')
!      .................... 0(n); 1 (p)
!       V_J is read from IMSRG and contains the following parts
!       nn/pp:  nlj1 <= nlj2
!       np:     No limitation on nlj1, nlj2
!      .............................
! ####  SINGLE-PARTICLE LABELS  ####
! ### index    n    l    2*j    2*tz      occ. 
!         1    0    0    1        -1    1
!      .............................
!      Read index for the s.p.s
!      .............................
       if(.NOT.allocated(sps_n)) ALLOCATE(sps_n(1:HO%nljmax*2))
       if(.NOT.allocated(sps_lj)) ALLOCATE(sps_lj(1:HO%nljmax*2))
       if(.NOT.allocated(sps_t)) ALLOCATE(sps_t(1:HO%nljmax*2))

          read(26,*) 
          read(26,*) 
       do tnlj=1,HO%nljmax*2   ! 2 is for two isospin projection
          read(26,*) id1,n1,l1,j1,tz1,occ 
          sps_n(id1) = n1  ! 0, 1, ... 
          sps_lj(id1)= (2*l1+j1-1)/2 
          sps_t(id1) = int(1-tz1)/2  ! transf: (-1 for p; +1 for n) to (1 for p; 0 for n) 
!          print *, id1,n1,l1,j1,tz1,sps_t(id1)
       enddo
          read(26,*) 
          read(26,*) 
          read(26,*) 
!####  TWO-BODY UNNORMALIZED MATRIX ELEMENTS  ####
!####  GAMMA((ab)J (cd)J) ==  GAMMA((cd)J (ab)J) 
!# a   b   c   d    2*J    GAMMA((ab)J (cd)J)
!4      read(26,*,end=5) t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,me
4      read(26,*,end=5) id1,id2,id3,id4,JJ,me  ! JJ was doubled in Nathan's code
       JJ   = JJ/2                           
       t1   = sps_t(id1)
       t2   = sps_t(id2)
       t3   = sps_t(id3)
       t4   = sps_t(id4)
       n1   = sps_n(id1)
       n2   = sps_n(id2)
       n3   = sps_n(id3)
       n4   = sps_n(id4)
       lj1  = sps_lj(id1)
       lj2  = sps_lj(id2)
       lj3  = sps_lj(id3)
       lj4  = sps_lj(id4)
       nlj1 = SPB%nlj(n1,lj1)   ! n1 =0, 1, 2  
       nlj2 = SPB%nlj(n2,lj2) 
       nlj3 = SPB%nlj(n3,lj3) 
       nlj4 = SPB%nlj(n4,lj4)

       if(it1234(t1,t2,t3,t4) .lt. 0 .or. it1234(t1,t2,t3,t4) .gt. 5) stop 'Error in it1234'

             l1   = SPB%l(nlj1) 
             l2   = SPB%l(nlj2) 
             l3   = SPB%l(nlj3) 
             l4   = SPB%l(nlj4) 

             twoj1 = SPB%twoj(nlj1) !  2*int(lj1/2) +1
             twoj2 = SPB%twoj(nlj2) ! 2*int(lj2/2) +1
             twoj3 = SPB%twoj(nlj3) ! 2*int(lj3/2) +1
             twoj4 = SPB%twoj(nlj4) ! 2*int(lj4/2) +1

       
             iphase_12 = iv((twoj1+twoj2)/2+JJ+1)
             iphase_34 = iv((twoj3+twoj4)/2+JJ+1)

         if(iCS.eq.1) then
            iphase_CS = 1
         else
            iphase_CS = iv((l1+l2-l3-l4)/2)
         endif

         V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4)) &
     &   = me*iphase_CS  !/H%ddd_tbme
!    ......... (12) -> (21)
         V_AB_CD_J(nlj2,nlj1,nlj3,nlj4,JJ,it1234(t2,t1,t3,t4)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12 
!    .........  (34) -> (43)
        V_AB_CD_J(nlj1,nlj2,nlj4,nlj3,JJ,it1234(t1,t2,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_34 
!    ......... (12) (34) -> (21) (43)
        V_AB_CD_J(nlj2,nlj1,nlj4,nlj3,JJ,it1234(t2,t1,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12*iphase_34 


!        write(300,'(4i3,4i4,5i5,f15.8)') t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,&
!      &        JJ,V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))
       goto 4
5      continue
!       write(*,*) ' The ME2B in J scheme is rewritten into',f2bJ_save
! ................................... transform from J to M 
       write(*,*) ' transform ME2B from J to M scheme'

      print *, ' uncoupling ME2B from J-scheme to m-scheme and store to',File%f2bm
        icheck = H%iabcd_max/10
        do iabcd = 1, H%iabcd_max
        if(mod(iabcd,icheck).eq.0) call progress(iabcd/icheck)
        ia = H%ka(iabcd)
        jja=tnljm%twoj(ia)  ! jang(ia)
        jma=tnljm%twom(ia)
        jta=tnljm%t(ia)      !mtisos(ia) 
        jna=tnljm%n(ia)
        jlja=tnljm%lj(ia)


        ib = H%kb(iabcd)
        jjb=tnljm%twoj(ib) ! jang(ib)
        jmb=tnljm%twom(ib)  ! mjang(ib)
        jtb=tnljm%t(ib)  !mtisos(ib)
        jnb=tnljm%n(ib)
        jljb=tnljm%lj(ib)

        delta_ab=0.d0
        if(jna.eq.jnb .and. jlja.eq.jljb) delta_ab=1.d0

        ic = H%kc(iabcd)
         jjc=tnljm%twoj(ic)   ! jang(ic)
         jmc=tnljm%twom(ic)   ! mjang(ic)
         jtc=tnljm%t(ic)      ! mtisos(ic)
         jnc=tnljm%n(ic)
         jljc=tnljm%lj(ic)
          id = H%kd(iabcd)
          jjd=tnljm%twoj(id)   !jang(id)
          jmd=tnljm%twom(id)   !mjang(id)
          jtd=tnljm%t(id)      !mtisos(id) -1 (p); +1 (n)
          jnd=tnljm%n(id)      ! 1,2,3
          jljd=tnljm%lj(id)

             
       nlj1 = SPB%nlj(jna-1,jlja)   ! n1 =0, 1, 2  
       nlj2 = SPB%nlj(jnb-1,jljb)
       nlj3 = SPB%nlj(jnc-1,jljc)
       nlj4 = SPB%nlj(jnd-1,jljd)
       t1   = mit(jta)  ! t1=0(n),1(p)
       t2   = mit(jtb)
       t3   = mit(jtc)
       t4   = mit(jtd)


!        if(jta+jtb.ne.0) cycle !!! jmyao: for test v_npnp


        delta_cd=0.d0
        if(jnc.eq.jnd .and. jljc.eq.jljd) delta_cd=1.d0

          suma=0.d0
          DO JJ=0,HO%twojmax  ! JJ is NOT doubled
!              if(JJ.ne.5) cycle !!! jmyao: for test
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ

              call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)
              call CJJ(jjc,jjd,2*JJ,jmc,jmd,2*MMJ,cb2)

!              cb1  = CG_Save((jja+1)/2,(jma+1)/2,(jjb+1)/2,(jmb+1)/2,JJ,MMJ) ! J12,M12 are CNOT doubled 
!              cb2  = CG_Save((jjc+1)/2,(jmc+1)/2,(jjd+1)/2,(jmd+1)/2,JJ,MMJ) ! J12,M12 are CNOT doubled 
!           ........................ normalization factor
!            phasJT=(-1.d0)**(JJ)
!            AN_AB=sqrt(1.d0-delta_ab*phasJT)/(1.d0+delta_ab)
!            AN_CD=sqrt(1.d0-delta_cd*phasJT)/(1.d0+delta_cd)
!
!            if(AN_AB.le.1e-15) cycle
!            if(AN_CD.le.1e-15) cycle
!
!            AN_INV=1.d0 !/(AN_AB*AN_CD)

!              cb3  = CG_Save(1,(jta+1)/2,1,(jtb+1)/2,IT,MMT) ! J12,M12 are NOT doubled 
!              cb4  = CG_Save(1,(jtc+1)/2,1,(jtd+1)/2,IT,MMT) ! J12,M12 are NOT doubled 

             sumando=cb1*cb2            &
      &      *V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))

             suma=suma+sumando

          END DO
          END DO
           H%ME2BM(iabcd) = suma
!           write(99) suma
       enddo ! iabcd 
!     ..............................................................................
      write(*,*) ' Writing ME2B in m-scheme into file: '
      do iabcd = 1, H%iabcd_max
         write(99) H%ME2BM(iabcd) ! suma
       enddo ! iabcd 
!  .........................
       end
