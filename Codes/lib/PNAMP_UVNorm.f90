      Subroutine normalization(ZV,ZU,ZNF,N)
!     ........................................
!     determination of normalization factor N
!     1=<a|a>=N^-1 * pf(...)
!     ........................................
      Use VAPHFB_Par
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)

      Dimension ZV(N,N),ZU(N,N)
      Dimension ZVC(N,N),ZUC(N,N)
      DIMENSION ZAUX1(N,N)
      DIMENSION ZAUX2(N,N)
      DIMENSION ZAUX3(N,N)
      DIMENSION ZAUX4(N,N)
      DIMENSION ZPfaf(N*2,N*2)
      DIMENSION Ipiv2(2*N,2)
      DIMENSION Ztest(2,2)
      DIMENSION Itest(4,2)
      NDIM = N

        do ii=1,NDIM
         do jj=1,NDIM
           ZUC(jj,ii)= dconjg(ZU(jj,ii))
           ZVC(jj,ii)= dconjg(ZV(jj,ii))
         end do
        end do

      ! zAUX1 = V^T * U  
        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZV,NDIM,ZU,NDIM,zzero,zAUX1,NDIM) 
      ! zAUX2 = V^T * V^*  
        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZV,NDIM,ZVC,NDIM,zzero,zAUX2,NDIM) 
      ! zAUX3 = -V^+ * V  
        call ZGEMM ('c','n',NDIM,NDIM,NDIM,zmone,ZV,NDIM,ZV,NDIM,zzero,zAUX3,NDIM) 
      ! zAUX4 = U^+ * V^* 
        call ZGEMM ('c','n',NDIM,NDIM,NDIM,zone,ZU,NDIM,ZVC,NDIM,zzero,zAUX4,NDIM) 

      !Big M
      do ii=1,N
       do jj=1,N
          ZPfaf(jj,ii)    = zAUX1(jj,ii)
          ZPfaf(jj+N,ii+N)= zAUX4(jj,ii)
!          ZPfaf(jj,ii+N)  = zAUX2(jj,ii)
!          ZPfaf(jj+N,ii)  = zAUX3(jj,ii)
           write(911,*) jj,ii,ZPfaf(jj,ii)
       end do
      end do
      call ZPfaffianF(ZPfaf,2*N,2*N,Ipiv2,ZdetPfaf)
      print *, ' norm=',ZdetPfaf
      ZNF = ZdetPfaf  
      ! ..... test
       Ztest(1,1) =  zzero
       Ztest(2,2) =  zzero
       Ztest(1,2) =  0.234 
       Ztest(2,1) = -0.234 
      call ZPfaffianF(Ztest,2,2,Itest,ZdetPfaf)
      print *, ' norm_test=',ZdetPfaf
 
      return
      end
