!          ..........................................................
!          <(12)J;TMT | V | (34)J;TMT>
!          ..........................................................
           subroutine GenerateV_JTMT()
           USE VAPHFB_PAR
           implicit none

           REAL(DP) CG
           integer t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ
           integer twoj1,twoj2,iphase_12   
           integer l1,l2,l3,l4
           integer twoj3,twoj4,iphase_34    
           real*8  fac_cg,me,cg_coeff(-1:1,-1:1,0:1),cb
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
!             call CJJ(1,1,2*TT,t1,t2, MMT,cb)
             cb = CG(1,t1,1,t2,2*TT,MMT)
             cg_coeff(t1,t2,TT) = cb
!             write(*,*)  t1,t2,TT,cb
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
       if(Input%IntType.eq.1) then
           File%f2bJ='../Int/IMSRG_'//Input%cFlow//'_'//trim(Input%cIntID)//&
        &  Input%ctpp//trim(Input%cValID)//'_'//Input%chwHO//'_J.dat'
       else
          File%f2bJ='../Int/IMSRG_'//Input%cFlow//'_'//trim(Input%cIntID)// &
        &  '_'//trim(Input%cValID)//'_J.dat'
       endif

        write(*,*) '    Read ME2B from ',File%f2bJ
        open(26,file=File%f2bJ,status='old')
!      .................... 0(n); 1 (p)
!       V_J is read from IMSRG and contains the following parts
!       nn/pp:  nlj1 <= nlj2
!       np:     No limitation on nlj1, nlj2
!      .............................
4      read(26,*,end=5) t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,me
       nlj1 = SPB%nlj(n1,lj1)   ! n1 =0, 1, 2  
       nlj2 = SPB%nlj(n2,lj2) 
       nlj3 = SPB%nlj(n3,lj3) 
       nlj4 = SPB%nlj(n4,lj4)
       write(66,'(9i4,f12.8)') t1,t2,t3,t4,nlj1,nlj2,nlj3,nlj4,JJ,me/H%ddd_tbme
       write(77,'(9i4,f12.8)') t1,t2,t3,t4,nlj1,nlj2,nlj3,nlj4,JJ,me

       if(it1234(t1,t2,t3,t4) .lt. 0 .or. it1234(t1,t2,t3,t4) .gt. 5) stop 'Error in it1234'

             l1   = SPB%l(nlj1) 
             l2   = SPB%l(nlj2) 
             l3   = SPB%l(nlj3) 
             l4   = SPB%l(nlj4) 

!             if(iv(l1+l2).ne.iv(l3+l4)) goto 4

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
     &   = me*iphase_CS/H%ddd_tbme
!    ......... (12) -> (21)
!         if(t1 .ne. t2 .or. nlj1 .ne. nlj2)  &
         V_AB_CD_J(nlj2,nlj1,nlj3,nlj4,JJ,it1234(t2,t1,t3,t4)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12 
!    .........  (34) -> (43)
!         if(t3 .ne. t4 .or. nlj3 .ne. nlj4)  &
        V_AB_CD_J(nlj1,nlj2,nlj4,nlj3,JJ,it1234(t1,t2,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_34 
!    ......... (12) (34) -> (21) (43)
!         if(t1 .ne. t2 .or. nlj1 .ne. nlj2 .or. t3 .ne. t4 .or. nlj3 .ne. nlj4)  &
        V_AB_CD_J(nlj2,nlj1,nlj4,nlj3,JJ,it1234(t2,t1,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12*iphase_34 

             if(JJ.eq.0 .and. nlj1.eq.1 .and.nlj2.eq.1 .and. nlj3.eq.2 .and. nlj4.eq.3 ) &
     &       write(*,'(5i3,f9.4)') t1,t2,t3,t4,JJ,V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))


       goto 4
5      continue
!       write(*,*) V_AB_CD_JTMT(1,1,2,3,6,1,-1)
! ................................... transform from J to JT
        write(*,*) ' transform V from J to JTT_z'

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

             if(nlj1.eq.1 .and.nlj2.eq.1 .and. nlj3.eq.2 .and. nlj4.eq.3 .and. JJ.eq.6 .and. TT.eq.1) &
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
       write(*,*) V_AB_CD_JTMT(1,1,2,3,6,1,-1)

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
