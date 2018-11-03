
       subroutine Read_Rho3B_J0(ie,Rho3BJJ,lpr)
       USE VAPHFB_Par
       implicit none
       logical lpr
       INTEGER ie,id123,id456
       real*8 cme
       real*8, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: Rho3BJJ  
!     ...............................
       print *, ' -> reading Rho3B ...'
10     read(ie,*,end=20) id123,id456,cme
          Rho3BJJ(id123,id456)=cme
        write(121,*) id123,id456,cme
       goto 10
20     continue
       if(lpr) call Check_Trace_Rho3B(Rho3BJJ,1)

       end


!       ...............................................
        subroutine Check_Trace_Rho3B(Rho3BJJ,IsRho)
        USE VAPHFB_PAR
        implicit none
        real*8 Trace(0:2),Rho3Bt
        real*8, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: Rho3BJJ  
        INTEGER id123,J123,t123,id456,nn,np,na,IsRho
!     ...............................
        if(IsRho.eq.1) &
     &  write(*,*) ' ... Check the Trace of Rho3B Matrix '
        if(IsRho.eq.0) &
     &  write(*,*) ' ... Check the Trace of L3B Matrix '


        Trace(0:2) = zero
        do id123=1,Rho3B%idx123_vmax
           J123 = Rho3B%J123(id123)
           t123 = Rho3B%t12(id123) + Rho3B%t3(id123)
!     ...................................... normalization
           Rho3Bt = Rho3BJJ(id123,id123)
           if(abs(Rho3Bt).le.1.d-5) cycle
           if(t123 .eq. 0 ) trace(0) = trace(0) + (J123+1.d0)*Rho3Bt
           if(t123 .eq. 3 ) trace(1) = trace(1) + (J123+1.d0)*Rho3Bt
           trace(2) = trace(2) + (J123+1.d0)*Rho3Bt
!           trace(2) = trace(2) + Rho3Bt
        if(IsRho.eq.0) write(900,'(2i4,2f12.8)') id123,t123, Rho3Bt,trace(2) 
        enddo

        if(IsRho.eq.1) then 
!       .......................
        nn = Nucl%nucleon(0) - Nucl%ncore(0)
        np = Nucl%nucleon(1) - Nucl%ncore(1)
        na = nn + np
        if(abs(Trace(0)-nn*(nn-1)*(nn-2)).lt.1.d-8) then
           write(*,*) ' Trace[Rho3B](nnn) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](nnn)=',Trace(0)
        endif
        if(abs(Trace(1)-np*(np-1)*(np-2)).lt.1.d-8) then
           write(*,*) ' Trace[Rho3B](ppp) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](ppp)=',Trace(1)
        endif

        if(abs(Trace(2)-na*(na-1)*(na-2)).lt.1.d-8) then
           write(*,*) ' Trace[Rho3B](ttt) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](ttt)=',Trace(2)
        endif
!       .......................
       else
         write(*,*) ' Trace[L3B](nnn)=',Trace(0)
         write(*,*) ' Trace[L3B](ppp)=',Trace(1)
         write(*,*) ' Trace[L3B](ttt)=',Trace(2)
       endif

       end


       subroutine Add_Rho3B_J0(iq1,iq2,Rho3BJJ,lpr)
        USE VAPHFB_PAR
        implicit none
        logical lpr
        integer   iq1,iq2
        real*8, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: Rho3BJJ  
        INTEGER id123,id456,n12
!     ...............................
!       n12 = 2/(1 + iq1/iq2)  ! 1 or 2
       n12 = 1 + iq1/iq2  ! 1 or 2
        do id123=1,Rho3B%idx123_vmax
        do id456=1,Rho3B%idx123_vmax
         GCMDens%Rho3BSum(id123,id456) = &
      &                             GCMDens%Rho3BSum(id123,id456) &
      &                           + GCM%FqJk(iq1,0,1)*GCM%FqJk(iq2,0,1)       &
      &                            *GCM%njkkqq(0,iq1,iq2)/n12         &
      &                            *(Rho3BJJ(id123,id456)+Rho3BJJ(id456,id123))

         enddo
         enddo

       end

        subroutine Check_Rho3B(lpr)
        USE VAPHFB_PAR
        implicit none
        logical    lpr
        if(lpr) then
          call Check_Trace_Rho3B(GCMDens%Rho3BSum,1)
          call Check_Symmetry_Rho3B(GCMDens%Rho3BSum,1)
        endif

        end

!    .................................................
        subroutine Rho3B2Lambda3B(ie,ZRho1BJ0,ZRho2BJJ,ZRho3BJJ,ZLambda3B)
!    .................................................
!    ZRho1BJ0 is defined in the whole space
!    ZRhO2BJJ is defined in the pf.val space
!    ZRho3BJJ is defined in the me3b.val space 
!    .................................................
!    ATTENTION: index for quantum numbers in diff. densities might be diff.
!    .................................................
        USE VAPHFB_PAR
        implicit none
        integer    ie
        real*8  znorm,Trace(0:1)
        real*8, Dimension(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3BJJ,ZLambda3B 
        real*8  ZRho2BJJ(-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        real*8  ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        real*8  Wigner_6j
        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER k1,t1,n1,lj1,l1,j1,twom1
        INTEGER k2,t2,n2,lj2,l2,j2,twom2
        INTEGER k3,t3,n3,lj3,l3,j3,twom3
        INTEGER k4,t4,n4,lj4,l4,j4,twom4
        INTEGER k5,t5,n5,lj5,l5,j5,twom5
        INTEGER k6,t6,n6,lj6,l6,j6,twom6
        INTEGER J123,P123,P456,J12,J45,P12,P45,t12,t45
        INTEGER id123,id456,b12,a12,b45,a45,M12,M45,M123
        INTEGER tnlj1,tnlj2,tnlj4,tnlj5,tnlj6
        INTEGER SPB_l_lj,SPB_twoj_lj
        INTEGER a4

        INTEGER iphase23,J23
        INTEGER fidx1,fidx2,fidx3,fidx4,fidx5,fidx6,fb23,fa23,fa45
        INTEGER iphase64,fa64,b64,P64,t64,J64,P23,t23
        INTEGER id1,id2,id3,id4,id5,id6
        INTEGER P56,t56,iphase56,J56,fa56
        INTEGER J31,fa31,fb45,iphase31,t31,P31,fb31
        INTEGER fb12,fa12,iphase12,iphase45
!       ........ start

        print *, ' ... Transform Rho3B to L3B '

        ZLambda3B = ZRho3BJJ
!      ............................ loop over id123 for lambda3B
        do id123=1,Rho3B%idx123_vmax
           J123 = Rho3B%J123(id123)            ! doubled

           J12 = Rho3B%J12(id123)
           P12 = Rho3B%P12(id123)
           t12 = Rho3B%t12(id123)
           t3  = Rho3B%t3(id123)
           lj3 = Rho3B%lj3(id123)
           n3  = Rho3B%n3(id123)
           l3  = SPB_l_lj(lj3)*2     ! doubled l3
           j3  = SPB_twoj_lj(lj3)    ! doubled
           P123 = mod(l3/2+P12,2)

           b12   = Rho3B%VTPB%block(J12,P12,t12)
           a12   = Rho3B%a12(id123)
!      ............. quantum numbers for 1 and 2
           n1  = Rho3B%VTPB%n1(b12,a12)
           t1  = Rho3B%VTPB%t1(b12,a12)
           L1  = Rho3B%VTPB%twol1(b12,a12)  ! doubled 
           j1  = Rho3B%VTPB%j1(b12,a12)  ! doubled
           lj1 = (L1+j1-1)/2
           n2  = Rho3B%VTPB%n2(b12,a12)
           t2  = Rho3B%VTPB%t2(b12,a12)
           L2  = Rho3B%VTPB%twol2(b12,a12)  ! doubled 
           j2  = Rho3B%VTPB%j2(b12,a12)  ! doubled
           lj2 = (L2+j2-1)/2
!      ............. check
           tnlj1 = Rho3B%VSPB%tnlj(lj1,n1,t1)
           tnlj2 = Rho3B%VSPB%tnlj(lj2,n2,t2)
!           if(a12 .ne. Rho3B%VTPB%a(b12,tnlj1,tnlj2)) &
!     &     stop 'Error in QNs for a(1,2)'
!      ............................ loop over id456 for lambda3B

        do id456=1,Rho3B%idx123_vmax

           if(J123.ne.Rho3B%J123(id456)) cycle

           J45 = Rho3B%J12(id456)
           P45 = Rho3B%P12(id456)
           t45 = Rho3B%t12(id456)
           t6  = Rho3B%t3(id456)
           lj6 = Rho3B%lj3(id456)
           n6  = Rho3B%n3(id456)
           l6  = SPB_l_lj(lj6)*2       ! doubled l6
           j6  = SPB_twoj_lj(lj6)
           P456 = mod(l6/2+P45,2)

           if(t12+t3.ne.t45+t6) cycle
           if(P123.ne.P456)     cycle

           b45   = Rho3B%VTPB%block(J45,P45,t45)
           a45   = Rho3B%a12(id456)
!      ............. quantum numbers for 1 and 2
           n4  = Rho3B%VTPB%n1(b45,a45)     ! n=0,1,2,3,...
           t4  = Rho3B%VTPB%t1(b45,a45)
           L4  = Rho3B%VTPB%twol1(b45,a45)  ! doubled 
           j4  = Rho3B%VTPB%j1(b45,a45)     ! doubled
           lj4 = (L4+j4-1)/2
           n5  = Rho3B%VTPB%n2(b45,a45)
           t5  = Rho3B%VTPB%t2(b45,a45)
           L5  = Rho3B%VTPB%twol2(b45,a45)  ! doubled 
           j5  = Rho3B%VTPB%j2(b45,a45)  ! doubled
           lj5 = (L5+j5-1)/2
           tnlj4 = Rho3B%VSPB%tnlj(lj4,n4,t4)
           tnlj5 = Rho3B%VSPB%tnlj(lj5,n5,t5)


!           if(id123.ne.4036 .or. id123.ne.4035) cycle
!           if(id123.eq.id456 .and. (id123.ge.4033 .and. id123.le.4038)) then 
!           write(*,*) 'id123=',id123,'J12=',J12,'J123=',J123
!           write(*,*) 't',t1,t2,t3 
!           write(*,*) 'n',n1,n2,n3 
!           write(*,*) 'l',l1/2,l2/2,l3/2
!           write(*,*) 'j',j1,j2,j3
!           endif
!       here

        t23 = t2 + t3
        P23 = mod((l2+l3)/2,2)

        t64 = t6+t4
        P64 = mod((l6+l4)/2,2)

        t31 = t1+t3
        P31 = mod((l1+l3)/2,2)

        t56 = t5+t6
        P56 = mod((l5+l6)/2,2)
        !if(id123.eq.2336 .and. id456.ne.2363) write(*,*)   
!  .... T1 
        if(lj1.eq.lj6 .and. t1.eq.t6 .and. P45.eq.P23 .and. t23.eq.t45) then

          J23 = J45
          if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else  ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
           if(t45.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb23 = TPB%block(J23,P23,t23)
           fa23 = TPB%a(fb23,fidx2,fidx3)
           fa45 = TPB%a(fb23,fidx4,fidx5)

           if(fb23.eq.-1 .or. fa23.eq.-1 .or. fa45.eq.-1) then
!           write(*,*) 'j',j2,j3,j4,j5
!           write(*,*) J23,P23,t23,fb23
!           write(*,*) fb23,fa23,fa64 
!           write(*,*) 't',t2,t3,t4,t5 
!           write(*,*) 'n',n2,n3,n4,n5 
!           write(*,*) 'l',l2,l3,l4,l5 
            ZRho2BJJ(fa23,fa45,fb23)=zero
           endif
!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J45+(j2+j3)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *iphase23*iphase45*ZRho2BJJ(fa23,fa45,fb23)      ! n=1,2,..

          endif

!  .... T2 

        if(lj1.eq.lj5 .and. t1.eq.t5 .and. P64.eq.P23 .and. t64.eq.t23) then

         do J23 = abs(j2-j3)/2, (j2+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(23)(64), if 2=3, or 4=6, J is only even
            if(t2.eq.t3 .and. lj2.eq.lj3 .and. iv(J23).eq.-1) cycle
            if(t4.eq.t6 .and. lj4.eq.lj6 .and. iv(J23).eq.-1) cycle
            if(J23.lt.abs(j6-j4)/2 .or. J23.gt.(j6+j4)/2) cycle

           J64 = J23
           if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else    ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t23.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t23.ne.1.and.id6.gt.id4) then
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j6+j4)/2+1)
            endif


           fb23 = TPB%block(J23,P23,t23)
           fa23 = TPB%a(fb23,fidx2,fidx3)
           fa64 = TPB%a(fb23,fidx6,fidx4)

           if(fa23.eq.-1 .or. fa64.eq.-1) then
!           write(*,'(10i6)') J23,P23,t23,fb23,fidx2,fidx3,fa23,fidx6,fidx4,fa64 
!           write(*,*) fb23,fa23,fa64 
!           write(*,*) 't',t2,t3,t6,t4 
!           write(*,*) 'n',n2,n3,n6,n4 
!           write(*,*) 'l',l2,l3,l6,l4 
           write(*,*) 'j',j2,j3,j6,j4
           stop
           endif

          ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv(J45+J23+(j1+j2+j3+j4)/2) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J23+1) &
     &                 *Wigner_6j(j4,j1,2*J45,J123,j6,2*J23) &
     &                 *Wigner_6j(j2,j3,2*J23,J123,j1,2*J12) &
     &                 *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                 *iphase23*iphase64*ZRho2BJJ(fa23,fa64,fb23)      ! n=1,2,..
          enddo
        endif

!  .... T3 

        if(lj1.eq.lj4 .and. t1.eq.t4 .and. &
     &     mod((l6+l5)/2,2) .eq. P23 .and. t23.eq.t5+t6) then

         do J23 = abs(j2-j3)/2, (j2+j3)/2 ! NOT Doubled

!        .......................................... Rho2B(23)(64), if 2=3, or 4=6, J is only even
            if(t2.eq.t3 .and. lj2.eq.lj3 .and. iv(J23).eq.-1) cycle
            if(t5.eq.t6 .and. lj5.eq.lj6 .and. iv(J23).eq.-1) cycle
            if(J23.lt.abs(j6-j5)/2 .or. J23.gt.(j6+j5)/2) cycle

           J56 = J23
           if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else    ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1
          endif

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t23.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t23.ne.1.and.id5.gt.id6) then
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j6+j5)/2+1)
            endif


           fb23 = TPB%block(J23,P23,t23)
           fa23 = TPB%a(fb23,fidx2,fidx3)
           fa56 = TPB%a(fb23,fidx5,fidx6)

           if(fa23.eq.-1 .or. fa56.eq.-1) then
           write(*,*) 'j',j2,j3,j5,j6
           stop
           endif
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j2+j3+j5+j6)/2) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J23+1) &
     &                 *Wigner_6j(j5,j1,2*J45,J123,j6,2*J23) &
     &                 *Wigner_6j(j2,j3,2*J23,J123,j1,2*J12) &
     &                 *ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                 *iphase23*iphase56*ZRho2BJJ(fa23,fa56,fb23)      ! n=1,2,..
          enddo
        endif

4     continue

!  .... T4 

        if(lj2.eq.lj6 .and. t2.eq.t6 .and. P45 .eq. P31 .and. t31.eq.t45) then

          J31 = J45
          if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
          else  ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id3.gt.id1) then
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb45 = TPB%block(J45,P45,t45)
           fa31 = TPB%a(fb45,fidx3,fidx1)
           fa45 = TPB%a(fb45,fidx4,fidx5)


           if(fb45.eq.-1 .or. fa31.eq.-1 .or. fa45.eq.-1) then
              ZRho2BJJ(fa31,fa45,fb45)=zero
           endif

          ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv((j1+j2)/2-J12+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *iphase31*iphase45*ZRho2BJJ(fa31,fa45,fb45)      ! n=1,2,..
          endif



!  .... T5 

        if(lj2.eq.lj5 .and. t2.eq.t5 .and. P64.eq.P31 .and. t31.eq.t64) then

         do J31 = abs(j1-j3)/2, (j1+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(31)(64), if 2=3, or 4=6, J is only even
            if(t1.eq.t3 .and. lj1.eq.lj3 .and. iv(J31).eq.-1) cycle
            if(t4.eq.t6 .and. lj4.eq.lj6 .and. iv(J31).eq.-1) cycle
            if(J31.lt.abs(j6-j4)/2 .or. J31.gt.(j6+j4)/2) cycle

           J64 = J31
           if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
         else    ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif

           if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t31.ne.1.and.id3.gt.id1) then
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t31.ne.1.and.id6.gt.id4) then
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j6+j4)/2+1)
            endif


           fb31 = TPB%block(J31,P31,t31)
           fa31 = TPB%a(fb31,fidx3,fidx1)
           fa64 = TPB%a(fb31,fidx6,fidx4)

           if(fa31.eq.-1 .or. fa64.eq.-1) then
           write(*,*) 'j',j3,j1,j6,j4
           stop
           endif
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j4+j1)/2+J12+J45+1) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J31+1) &
     &                 *Wigner_6j(j4,j2,2*J45,J123,j6,2*J31) &
     &                 *Wigner_6j(j1,j3,2*J31,J123,j2,2*J12) &
     &                 *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                 *iphase31*iphase64*ZRho2BJJ(fa31,fa64,fb31)      ! n=1,2,..
          enddo
        endif


!  .... T6 
        if(lj2.eq.lj4 .and. t2.eq.t4 .and. P56 .eq. P31 .and. t31.eq.t56) then

         do J31 = abs(j1-j3)/2, (j1+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(31)(64), if 2=3, or 4=6, J is only even
            if(t1.eq.t3 .and. lj1.eq.lj3 .and. iv(J31).eq.-1) cycle
            if(t5.eq.t6 .and. lj5.eq.lj6 .and. iv(J31).eq.-1) cycle
            if(J31.lt.abs(j6-j5)/2 .or. J31.gt.(j6+j5)/2) cycle

           J56 = J31
           if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
          else    ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif

           if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1
          endif

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t31.ne.1.and.id3.gt.id1) then
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t31.ne.1.and.id5.gt.id6) then
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j6+j5)/2+1)
            endif


           fb31 = TPB%block(J31,P31,t31)
           fa31 = TPB%a(fb31,fidx3,fidx1)
           fa56 = TPB%a(fb31,fidx5,fidx6)

           if(fa31.eq.-1 .or. fa56.eq.-1) then
           write(*,*) 'j',j3,j1,j5,j6
           stop
           endif
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j1+j2+j5+j6)/2+J12+J31) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J31+1) &
     &                 *Wigner_6j(j5,j2,2*J45,J123,j6,2*J31) &
     &                 *Wigner_6j(j1,j3,2*J31,J123,j2,2*J12) &
     &                 *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                 *iphase31*iphase56*ZRho2BJJ(fa31,fa56,fb31)      ! n=1,2,..
          enddo
        endif



!  .... T7 

        if(lj3.eq.lj6 .and. t3.eq.t6 .and. J12.eq.J45 .and. P12.eq.P45 .and. t12.eq.t45) then


!         .............(1,2) only contains nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

!         .............(4,5) only contains nn,np,pp
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
!      ................(1,2)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id1.gt.id2) then
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j1+j2)/2+1)
            endif
!      ................(4,5)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb45 = TPB%block(J45,P45,t45)
           fa12 = TPB%a(fb45,fidx1,fidx2)
           fa45 = TPB%a(fb45,fidx4,fidx5)

           if(fb45.eq.-1 .or. fa12.eq.-1 .or. fa45.eq.-1) then
              ZRho2BJJ(fa12,fa45,fb45)=zero
           endif

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - ZRho1BJ0(t3,lj3,n3,n6)               & ! n=0,1,2,..
     &                          *iphase12*iphase45*ZRho2BJJ(fa12,fa45,fb45)      ! n=1,2,..
          endif


8      continue
!  .... T8 
        if(lj3.eq.lj5 .and. t3.eq.t5 .and. P12.eq.P64  .and. t12.eq.t64) then

          J64 = J12
          if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id1.gt.id2) then
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j2+j1)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id6.gt.id4) then
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j4+j6)/2+1)
            endif

!           do J23=abs(j2-j3)/2, (j2+j3)/2   ! NOT doubled 
           fb12 = TPB%block(J12,P12,t12)
           fa12 = TPB%a(fb12,fidx1,fidx2)
           fa64 = TPB%a(fb12,fidx6,fidx4)

           if(fb12.eq.-1 .or. fa12.eq.-1 .or. fa64.eq.-1) then
              ZRho2BJJ(fa12,fa64,fb12)=zero
           endif

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J45+(j3+j4)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j4,j3,2*J45,J123,j6,2*J12) &
     &                          *ZRho1BJ0(t3,lj3,n3,n5)               & ! n=0,1,2,..
     &                          *iphase12*iphase64*ZRho2BJJ(fa12,fa64,fb12)      ! n=1,2,..
          endif


9      continue


!  .... T9 
        if(lj3.eq.lj4 .and. t3.eq.t4 .and. P12.eq.P56  .and. t12.eq.t56) then

          J56 = J12
          if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1

         endif

           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id1.gt.id2) then
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j2+j1)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id5.gt.id6) then
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j5+j6)/2+1)
            endif
           fb12 = TPB%block(J12,P12,t12)
           fa12 = TPB%a(fb12,fidx1,fidx2)
           fa56 = TPB%a(fb12,fidx5,fidx6)


           if(fb12.eq.-1 .or. fa12.eq.-1 .or. fa56.eq.-1) then
              ZRho2BJJ(fa12,fa56,fb12)=zero
           endif
!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J12+(j5+j6)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j5,j3,2*J45,J123,j6,2*J12) &
     &                          *ZRho1BJ0(t3,lj3,n3,n4)               & ! n=0,1,2,..
     &                          *iphase12*iphase56*ZRho2BJJ(fa12,fa56,fb12)      ! n=1,2,..
          endif



10      continue
!  .... T10 
        if( lj1.eq.lj6 .and. t1.eq.t6 .and. lj2.eq.lj5 .and. t2.eq.t5 .and. lj3.eq.lj4 .and. t3.eq.t4) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         +2*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n4)                 ! n=0,1,2,..
          endif
!  .... T11 
        if( lj1.eq.lj6 .and. t1.eq.t6 .and. lj2.eq.lj4 .and. t2.eq.t4 .and. lj3.eq.lj5 .and. t3.eq.t5) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         +2*iv(J45+(j2+j3)/2+1)                 &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n5)                 ! n=0,1,2,..
          endif


!  .... T12 
        if( lj1.eq.lj5 .and. t1.eq.t5 .and. lj2.eq.lj6 .and. t2.eq.t6 .and. lj3.eq.lj4 .and. t3.eq.t4) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         +2*iv((j1+j2)/2-J12+1)                 &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n4)                 ! n=0,1,2,..
          endif
!  .... T13 
        if( lj1.eq.lj5 .and. t1.eq.t5 .and. lj2.eq.lj4 .and. t2.eq.t4 .and. lj3.eq.lj6 .and. t3.eq.t6) then

        if(J12.eq.J45) ZLambda3B(id123,id456) = ZLambda3B(id123,id456)&
     &                         -2*iv((j1+j2)/2-J12)                   &
     &                          *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n6)                 ! n=0,1,2,..
          endif

14      continue
!  .... T14 
        if( lj1.eq.lj4 .and. t1.eq.t4 .and. lj2.eq.lj6 .and. t2.eq.t6 .and. lj3.eq.lj5 .and. t3.eq.t5) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         -2*iv((j3+j2)/2-J12+J45)               &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n5)                 ! n=0,1,2,..
          endif

15      continue
!  .... T15 
        if( lj1.eq.lj4 .and. t1.eq.t4 .and. lj2.eq.lj5 .and. t2.eq.t5 .and. lj3.eq.lj6 .and. t3.eq.t6) then

        if(J12.eq.J45) ZLambda3B(id123,id456) = ZLambda3B(id123,id456)&
     &                        +2*ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n6)                 ! n=0,1,2,..
          endif
!  ..............................................
           write(333,'(2i6,f20.10)') id123,id456,ZLambda3B(id123,id456)
           enddo
           enddo
        end
!   ........................................
        subroutine Write_L3B(ie,Lambda3B,lpr)
        USE VAPHFB_PAR
        implicit none
        logical lpr
        integer ie
        real*8 L3B
        real*8 Lambda3B(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) 
        INTEGER id123,id456
!     ...............................
        open(ie,file='GCM_L3BSum.dat',status='unknown')
        do id123=1,Rho3B%idx123_vmax
        do id456=1,Rho3B%idx123_vmax

           L3B = Lambda3B(id123,id456)
!     ...................................... normalization
!           if(t123.eq.0 .and. id123.eq.id456) trace(0) = trace(0) + (J123+1.d0)*L3B
!           if(t123.eq.3 .and. id123.eq.id456) then
!            trace(1) = trace(1) + (J123+1.d0)*L3B
!           endif

!           if(id123.eq.id456) &
!        &  write(911,'(i8,f12.8)') id123, L3B
!       &      trace(2) = trace(2) + (J123+1.d0)*L3B 

           !if(abs(L3B).gt.CHOP) write(ie,'(2i8,f20.10)') id123,id456, L3B
           if(abs(L3B).gt.1.d-5) write(ie,'(2i8,f20.10)') id123,id456, L3B
        enddo
        enddo
       if(lpr) call Check_Trace_Rho3B(GCMDens%L3BSum,0)
       if(lpr) call Check_Symmetry_Rho3B(GCMDens%L3BSum,0)
       end


!      ...........................................
        subroutine Check_Symmetry_Rho3B(R3B,IsRho)
!      ...........................................
        USE VAPHFB_PAR
        implicit none
        integer IsRho
        real*8 R3B(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax)
        real*8 diff
        INTEGER id1m,id2m,id123,id456
!     ...............................
        if(IsRho.eq.1) &
     &  write(*,*) ' ... Check the Symmetry of Rho3B Matrix  '
        if(IsRho.eq.0) &
     &  write(*,*) ' ... Check the Symmetry of L3B Matrix  '
        diff=zero
        id1m = 0 
        id2m = 0 
!     ...............................
        do id123=1,Rho3B%idx123_vmax
        do id456=1,Rho3B%idx123_vmax

           if(abs(R3B(id123,id456) - R3B(id456,id123)).gt.abs(diff)) then 
             diff = R3B(id123,id456) - R3B(id456,id123)
             id1m = id123
             id2m = id456
           endif
        enddo
        enddo
        call print_basis(id1m,id2m)
        write(*,*) 'Symmetry is satisfied within the precision:',diff
        write(*,*) 'The largest difference is found for',id1m,id2m
        end


        subroutine print_basis(id123,id456)
        USE VAPHFB_PAR
        implicit none
        integer id123,id456
        integer J123,P123,b12,a12,J12,P12,t12
        integer t1,lj1,n1,l1,j1
        integer t2,lj2,n2,l2,j2
        integer t3,lj3,n3,l3,j3

        integer b45,a45,J45,P45,t45
        integer t4,lj4,n4,l4,j4
        integer t5,lj5,n5,l5,j5
        integer t6,lj6,n6,l6,j6
        INTEGER SPB_l_lj,SPB_twoj_lj
           J123 = Rho3B%J123(id123)            ! doubled

           J12 = Rho3B%J12(id123)
           P12 = Rho3B%P12(id123)
           t12 = Rho3B%t12(id123)
           t3  = Rho3B%t3(id123)
           lj3 = Rho3B%lj3(id123)
           n3  = Rho3B%n3(id123)
           l3  = SPB_l_lj(lj3)*2     ! doubled l3
           j3  = SPB_twoj_lj(lj3)    ! doubled
           P123 = mod(l3/2+P12,2)

           b12   = Rho3B%VTPB%block(J12,P12,t12)
           a12   = Rho3B%a12(id123)
!      ............. quantum numbers for 1 and 2
           n1  = Rho3B%VTPB%n1(b12,a12)
           t1  = Rho3B%VTPB%t1(b12,a12)
           L1  = Rho3B%VTPB%twol1(b12,a12)  ! doubled 
           j1  = Rho3B%VTPB%j1(b12,a12)  ! doubled
           lj1 = (L1+j1-1)/2
           n2  = Rho3B%VTPB%n2(b12,a12)
           t2  = Rho3B%VTPB%t2(b12,a12)
           L2  = Rho3B%VTPB%twol2(b12,a12)  ! doubled 
           j2  = Rho3B%VTPB%j2(b12,a12)  ! doubled
           lj2 = (L2+j2-1)/2

           J45 = Rho3B%J12(id456)
           P45 = Rho3B%P12(id456)
           t45 = Rho3B%t12(id456)
           t6  = Rho3B%t3(id456)
           lj6 = Rho3B%lj3(id456)
           n6  = Rho3B%n3(id456)
           l6  = SPB_l_lj(lj6)*2       ! doubled l6
           j6  = SPB_twoj_lj(lj6)

           b45   = Rho3B%VTPB%block(J45,P45,t45)
           a45   = Rho3B%a12(id456)
!      ............. quantum numbers for 1 and 2
           n4  = Rho3B%VTPB%n1(b45,a45)     ! n=0,1,2,3,...
           t4  = Rho3B%VTPB%t1(b45,a45)
           L4  = Rho3B%VTPB%twol1(b45,a45)  ! doubled 
           j4  = Rho3B%VTPB%j1(b45,a45)     ! doubled
           lj4 = (L4+j4-1)/2
           n5  = Rho3B%VTPB%n2(b45,a45)
           t5  = Rho3B%VTPB%t2(b45,a45)
           L5  = Rho3B%VTPB%twol2(b45,a45)  ! doubled 
           j5  = Rho3B%VTPB%j2(b45,a45)  ! doubled
           lj5 = (L5+j5-1)/2
           write(*,*) 'J123=',J123,' J12=',J12,' J45=',J45
           write(*,10) '(123):', n1,l1,j1,n2,l2,j2,n3,l3,j3
           write(*,10) '(456):', n4,l4,j4,n5,l5,j5,n6,l6,j6
 10   format(a,3i4,',',3i4,',',3i4)
           end
