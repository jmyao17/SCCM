!          ..........................................................
!          read ME2B in J scheme from the file 'File%f2bJ'
!          store ME2B in m scheme into the file 'File%f2bm'
!          ..........................................................
!          <(12)J | V | (34)J> ----> <12| V |34>
!          ..........................................................
!          JJmax - > J12max
!          ..........................................................
           subroutine INT_ME2B_J2M()
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
           integer J12max,t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ
           integer twoj1,twoj2,iphase_12   
           integer l1,l2,l3,l4
           integer twoj3,twoj4,iphase_34    
           real*8  fac_cg,me,cg_coeff(-1:1,-1:1,0:1),cb
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


        open(26,file=trim(INT_DIR)//trim(File%f2bJ),status='old')
        open(99,file=trim(INT_DIR)//trim(File%f2bm),form='unformatted',status='unknown')

!     .............. allocation and initialization

       if(.NOT. ALLOCATED(V_AB_CD_J))  &
     & ALLOCATE(V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5))
       V_AB_CD_J(HO%nljmax,HO%nljmax,HO%nljmax,HO%nljmax,0:HO%TwoJmax,0:5) = zero 

!      ............. proton(-1); neutron (+1)
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
         ! initialization
         it1234(it1,it2,it3,it4) = -1
         if(it1+it2 .ne. it3+it4) cycle
         it1234(it1,it2,it3,it4) = ii  ! 0 - 5
         itt1(ii) = it1
         itt2(ii) = it2
         itt3(ii) = it3
         itt4(ii) = it4
         ii = ii + 1
      enddo
      enddo
      enddo
      enddo
!      only two-body
       if(Input%IntType.eq.1) then
          if(Input%InME3B.ne.1) & 
        &  File%f2bJ='../Int/IMSRG_'//trim(Input%cFlow)//'_'//Input%cIntID//&
        &  Input%ctpp//Input%cValID//'_'//Input%chwHO//'_J.dat'
          if(Input%InME3B.eq.1) & 
        &  File%f2bJ='../Int/IMSRG_'//trim(Input%cFlow) &
        &  //'_chi2b00000_srg0625'//&
        &  Input%ctpp//Input%cValID//'_'//Input%chwHO//'_J.dat'
       else
          File%f2bJ='../Int/IMSRG_'//trim(Input%cFlow)//'_'&
     &    //trim(Input%cIntID)//'_'//trim(Input%cValID)//'_J.dat'
       endif


        call StripSpaces(File%f2bJ)
        write(*,'(a20,a80)') '    Read ME2B from ',trim(File%f2bJ)
!      .................... 0(n); 1 (p)
!       V_J is read from IMSRG and contains the following parts
!       nn/pp:  nlj1 <= nlj2
!       np:     No limitation on nlj1, nlj2
!      .............................
       J12max = 0



4      read(26,*,end=5) t1,t2,t3,t4,n1,n2,n3,n4,lj1,lj2,lj3,lj4,JJ,me
       nlj1 = SPB%nlj(n1,lj1)   ! n1 =0, 1, 2  
       nlj2 = SPB%nlj(n2,lj2) 
       nlj3 = SPB%nlj(n3,lj3) 
       nlj4 = SPB%nlj(n4,lj4)
       if(JJ .gt. J12max) J12max=JJ

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
     &   = me*iphase_CS !/H%ddd_tbme already taken into account in me

!    ......... (12)(34) -> (34)(12): included in the TBME files
         V_AB_CD_J(nlj3,nlj4,nlj1,nlj2,JJ,it1234(t3,t4,t1,t2)) &
     &   =V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4)) 


!    ......... (12) -> (21)
        V_AB_CD_J(nlj2,nlj1,nlj3,nlj4,JJ,it1234(t2,t1,t3,t4)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12 
!    .........  (34) -> (43)
        V_AB_CD_J(nlj1,nlj2,nlj4,nlj3,JJ,it1234(t1,t2,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_34 
!    ......... (12) (34) -> (21) (43)
        V_AB_CD_J(nlj2,nlj1,nlj4,nlj3,JJ,it1234(t2,t1,t4,t3)) &
     & = V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))*iphase_12*iphase_34 



       goto 4
5      continue
! ................................... transform from J to M
      write(*,'(a30)') '     Transform NN from J to M: ' 
!     only decouple the nonzero matrix elements 
        icheck = H%iabcd_max/100

      do iabcd = 1, H%iabcd_max
         if(mod(iabcd,icheck).eq.0) call progress_bar(iabcd/icheck) ! call progress(iabcd/icheck)
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


          suma=0.d0
          DO JJ=0,J12max !HO%twojmax  ! JJ is NOT doubled
           DO MJ=0,2*JJ
              MMJ=-JJ+MJ
              cb1  = CG_Save((jja+1)/2,(jma+1)/2,(jjb+1)/2,(jmb+1)/2,JJ,MMJ) ! J12,M12 are CNOT doubled 
              cb2  = CG_Save((jjc+1)/2,(jmc+1)/2,(jjd+1)/2,(jmd+1)/2,JJ,MMJ) ! J12,M12 are CNOT doubled 

             sumando=cb1*cb2            &
      &      *V_AB_CD_J(nlj1,nlj2,nlj3,nlj4,JJ,it1234(t1,t2,t3,t4))

             suma=suma+sumando

          END DO
          END DO
           H%ME2BM(iabcd) = suma
!           write(99) suma
       enddo ! iabcd 
!     ..............................................................................
      DEALLOCATE(V_AB_CD_J)  
!     ..............................................................................
      write(*,*)
      write(*,'(a30)',advance='no') &
           '     Writing ME2B in m-scheme into file: '
      icheck = H%iabcd_max/100
      do iabcd = 1, H%iabcd_max
         if(mod(iabcd,icheck).eq.0) call progress_bar(int(iabcd/icheck))
         write(99) H%ME2BM(iabcd) ! suma
       enddo ! iabcd    
      write(*,*)
 
      end
