!          ................................................................
!          Read the ME2B files (in J-scheme) from Nathan and transform into
!          <(12)J;TMT | V | (34)J;TMT>
!          .................................................................
           subroutine Generate_ME2B()
           USE VAPHFB_PAR
           implicit none

           REAL(DP) CG
           integer tnlj,tz1
           integer t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ
           integer twoj1,twoj2,iphase_12   
           integer l1,l2,l3,l4
           integer j1,id1,id2,id3,id4
           integer twoj3,twoj4,iphase_34    
           real*8  fac_cg,me,cg_coeff(-1:1,-1:1,0:1),cb,occ
           real*8  VJTA
           real*8, dimension(:,:,:,:,:,:,:), allocatable :: V_AB_CD_JTMT
           real*8, dimension(:,:,:,:,:,:), allocatable :: V_AB_CD_J
           integer ia,ib,ic,id,iT,iJ,iphase_cd
           integer iphaseJ,iphaseJT
           integer nlj1,nlj2,nlj3,nlj4,MMT,MT,TT
           real*8  delta_J, delta_JT,d_ab,d_cd,d_abcd
           real*8  d_12,d_34,d_1234
           real*8  AN_12,AN_34,AN_AB,AN_CD
           integer it1234(0:1,0:1,0:1,0:1)
           integer itt1(0:5),itt2(0:5),itt3(0:5),itt4(0:5)
           integer ii,it1,it2,it3,it4
           integer iphase_CS
           integer iprint
           real*8 diff
           character(len=100) :: f2bJ_save

           integer, dimension(:), allocatable :: sps_t,sps_n,sps_lj

!     .............. initialization
       if(.NOT. ALLOCATED(V_AB_CD_J))  &
     & ALLOCATE(V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5))
       if(.NOT. ALLOCATED(V_AB_CD_JTMT))  &
     & ALLOCATE(V_AB_CD_JTMT(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:1,-1:1))
       do iT=0,HO%Tmax
       do MMT=-iT,iT
       do ia=1,HO%nljmax
       do ib=1,HO%nljmax
        do ic=1,HO%nljmax
         do id=1,HO%nljmax
          do iJ=0,HO%TwoJmax
!  ...................
                 V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT)=0.d0
!  ...................
           end do
           end do
          end do
         end do
        end do
       end do
      end do

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

          File%f2bJ=File%IMSRG_Hme2b 

!          if(Input%InME3B.eq.1) & 
!        &  File%f2bJ='../Int/IMSRG_'//Input%cFlow//'_chi2b00000_srg0625'//&
!        &  Input%ctpp//Input%cValID//'_'//Input%chwHO//'_J.dat'
!       else
          f2bJ_save='../Int/IMSRG_'//Input%cFlow//'_'&
     &    //Input%cIntID//'_J.dat'
!       endif

        open(300,file=f2bJ_save,status='unknown')

        write(*,*) '    Read ME2B from ',File%f2bJ
        open(26,file=File%f2bJ,status='old')
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


        write(300,'(4i3,4i4,5i5,f15.8)') t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,&
      &        JJ,V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))
       goto 4
5      continue
       write(*,*) ' The ME2B in J scheme is rewritten into',f2bJ_save
! ................................... transform from J to JT
        write(*,*) ' transform ME2B from J to JTT_z'

       do nlj1=1,HO%nljmax
       do nlj2=1,HO%nljmax

       do nlj3=1,HO%nljmax
       do nlj4=1,HO%nljmax

!          if(nlj3+nlj4 .lt. nlj1+nlj2) cycle

          do JJ=0,HO%TwoJmax
            do TT=0,HO%Tmax
            do ii=0,5
            t1 = itt1(ii)  ! 0 (n); 1 (p)
            t2 = itt2(ii)
            t3 = itt3(ii)
            t4 = itt4(ii)

            MMT = (iv(t1)+iv(t2))/2   ! iv(0)= 1 for n; iv(1)=-1 for p;
            if(t1+t2 .ne. t3+t4) cycle

            delta_J =one ! /(2*TT+1) ! JTM_T (degeneracy: 2*TT+1)
! .............................................
!     ........ for nnnn and pppp        
                  V_AB_CD_JTMT(nlj1,nlj2,nlj3,nlj4,JJ,TT,MMT)=  &
     &            V_AB_CD_JTMT(nlj1,nlj2,nlj3,nlj4,JJ,TT,MMT)   & 
     &            +V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))  &
     &            *delta_J*cg_coeff(iv(t1),iv(t2),TT)*cg_coeff(iv(t3),iv(t4),TT)

             if(nlj1.eq.1 .and.nlj2.eq.1 .and. nlj3.eq.1 .and. nlj4.eq.1 .and. JJ.eq.1 .and. TT.eq.0) &
     &       write(*,'(6i3,4f9.4)') t1,t2,t3,t4,TT,MMT,V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4)),&
     &                               V_AB_CD_JTMT(nlj1,nlj2,nlj3,nlj4,JJ,TT,MMT),   &
     &             cg_coeff(iv(t1),iv(t2),TT)*cg_coeff(iv(t3),iv(t4),TT),delta_J                      
              enddo ! ii
! ...................................
        enddo
        enddo
! ...............
        enddo
        enddo
        enddo
        enddo
!      ................ print
        write(*,*) '    Print ME2B to ',File%f2bJTMT
        open(36,file=File%f2bJTMT,status='unknown')
       do ia=1,HO%nljmax
       do ib=ia,HO%nljmax
                     d_ab = zero
        if(ia.eq.ib) d_ab = one 
        do ic=1,HO%nljmax
         do id=ic,HO%nljmax

            if(ic+id .lt. ia+ib) cycle
            if(ic+id .eq. ia+ib .and. ic.lt.ia) cycle

                           d_cd = zero
              if(ic.eq.id) d_cd = one
!              d_abcd =  (1.d0+d_ab)*(1.d0+d_cd) 
               
          do iJ=0,HO%TwoJmax
          do iT=0,HO%Tmax

            iphaseJT  = iv(iJ+iT)
!            iphase_cd = iv((twoj3+twoj4)/2+iJ+iT)
            AN_AB     = sqrt((1.d0-d_ab*iphaseJT))/(1.d0+d_ab)
            AN_CD     = sqrt((1.d0-d_cd*iphaseJT))/(1.d0+d_cd)
            delta_JT  = AN_AB*AN_CD
            if(abs(delta_JT).le.CHOP) cycle

            do MMT=-iT,iT
               VJTA   = V_AB_CD_JTMT(ia,ib,ic,id,iJ,iT,MMT)*delta_JT   
!               V_AB_CD_JT(ia,ib,ic,id,iJ,iT) = VJTA

!            if(ia.eq.1 .and. ib.eq.2 .and. ic.eq.1 .and. id.eq.4)
            if(abs(VJTA) .gt. CHOP)  &
     &      write(36,'(7i5,f15.8)') ia,ib,ic,id,iJ,iT,MMT,VJTA
!             if(ia.eq.1 .and.ib.eq.4 .and. ic.eq.2 .and. id.eq.2 .and. JJ.eq.3 .and. TT.eq.0) &
           end do  ! MMT

           end do
          end do
         end do
        end do
       end do
      end do
!  .........................
       end
