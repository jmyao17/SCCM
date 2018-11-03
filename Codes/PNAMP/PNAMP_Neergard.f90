      Subroutine Neegard_Wust(ZU0,ZV0,ZU1bar,ZV1bar,ZZTZ1,N,znorm)
      Use VAPHFB_Par
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      integer ii,jj

      Dimension ZV0(N,N),ZU0(N,N)
      Dimension ZV1bar(N,N),ZU1bar(N,N)
      Dimension ZU0INV(N,N),ZZ0(N,N),ZZ0c(N,N)
      Dimension ZU1barINV(N,N),ZZ1(N,N),ZZ1c(N,N)
      Dimension ZZTZ1(N,N)
      Dimension IPIV(N),IPIV1(N)
      Dimension ZWORK(N),ZWORK1bar(N)
      do ii=1,N
       do jj=1,N
        ZU0INV(jj,ii)   =ZU0(jj,ii)
        ZU1barINV(jj,ii)=ZU1bar(jj,ii)
       end do
      end do
!      print *, 'computing U^-1_0'
!               ================================
! ................ compute inverse of U0
!               ================================
      call ZGETRF(N,N,ZU0INV,N,IPIV,INFO1)
      call ZGETRI(N,ZU0INV,N,IPIV,ZWORK,N,INFO2)
      call ZGEMM ('n','n',N,N,N,zone,ZV0,N,ZU0INV,N,zzero,ZZ0,N)
      do ii=1,N
      do jj=1,N
         ZZ0c(jj,ii)=dconjg(ZZ0(jj,ii))
      end do
      end do
!               ================================
!   ............ compute the inverse of U1bar
!               ================================
      call ZGETRF(N,N,ZU1barINV,N,IPIV1,INFO1)
      call ZGETRI(N,ZU1barINV,N,IPIV1,ZWORK1bar,N,INFO2)
      call ZGEMM ('n','n',N,N,N,zone,ZV1bar,N,ZU1barINV,N,zzero,ZZ1,N)

      do ii=1,N
      do jj=1,N
         ZZ1c(jj,ii)=dconjg(ZZ1(jj,ii))
      end do
      end do
! ................ compute M=Z^t_1 * Z_0 ?= Z_0 * Z^t_1
!      print *, '  computing M '
      call ZGEMM ('t','n',N,N,N,zone,ZZ0c,N,ZZ1c,N,zzero,ZZTZ1,N)
!      call ZGEMM ('n','c',N,N,N,zone,ZZ1c,N,ZZ0c,N,zzero,ZZTZ1,N)
      return
      end

      subroutine zdiag(N,ZZTZ1,znorm)
      Use VAPHFB_Par
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)

      Dimension ZZTZ1(N,N)
      Dimension ZWR(N), ZWR2(N)
      Dimension ZVL(N,N)
      Dimension ZVR(N,N)
      Dimension ZWORK1(2*N)
      real*8 WORK1(2*N)

!      ZZTZ1 = zzero
!      open(800,file='fort.800',status='old')
      do jj=1,N
      do ii=1,N
        write(800,*) ii,jj,ZZTZ1(ii,jj)
       enddo 
       enddo 
      call ZGEEV('n','n',N,ZZTZ1,N,ZWR,ZVL,N,ZVR,N,ZWORK1,2*N,WORK1,INFO3)
      ztest = zzero
      zdet  = zone
      znorm = zzero
!    ........................ ordering the eigen values, which are two-fold degenerate
      call zordls(N,ZWR)

!     judge if two- or four-fold degenarate
      ii = 1   
      do while(ii.lt.N)
          if(abs(dreal(ZWR(ii))-dreal(ZWR(ii+1)))/abs(dreal(ZWR(ii))) .lt. 0.01d0 .and. &
          &  abs(dreal(ZWR(ii))-dreal(ZWR(ii+2)))/abs(dreal(ZWR(ii))) .ge. 0.01d0 ) then
            ii         = ii+2
          endif
!     ......................
          if(abs(dreal(ZWR(ii))-dreal(ZWR(ii+1)))/abs(dreal(ZWR(ii))) .lt. 0.01d0  .and. &
           & abs(dreal(ZWR(ii))-dreal(ZWR(ii+2)))/abs(dreal(ZWR(ii))) .lt. 0.01d0 ) then
!     ......................
            do ll=ii,ii+3
               kk         = ll
               zp         = ZWR(ll)
               if(ll.eq.ii+3) cycle 
               do 20 jj = ll+1,ii+3
                 if (dimag(ZWR(jj)).gt.dimag(zp)) then
                    kk = jj
                    zp = ZWR(jj)
                endif
   20          continue
               if (kk.ne.ll) then  ! exchange k, l
                  ZWR (kk)  = ZWR(ll)
                  ZWR(ll)  = zp
               endif                                     
             enddo
             ii     = ii+4
          endif
       enddo

       do ii=1,N,2 !2
!           print *,ZWR(ii)
          ztest = ztest + cdlog(zone+ZWR(ii))  ! log(x): the logarithm to the base e.
       enddo
       znorm  = CDEXP(ztest)
!       znorm  = CDEXP(ztest/2)
      return
      end

