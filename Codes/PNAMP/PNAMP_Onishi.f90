      subroutine Onishi(ZU0,ZV0,ZU1bar,ZV1bar,znorm,NDIM)           
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))

	DIMENSION ZU0(NDIM,NDIM),ZV0(NDIM,NDIM)
	DIMENSION ZU1bar(NDIM,NDIM),ZV1bar(NDIM,NDIM)
	DIMENSION ZAUX1(NDIM,NDIM)
	DIMENSION ZAUX2(NDIM,NDIM),ZX(NDIM,NDIM)
	DIMENSION ZAUX3(NDIM,NDIM)
	DIMENSION ZAUX4(NDIM,NDIM)	
	DIMENSION ZAUX5(NDIM,NDIM)
	DIMENSION ZAUX6(NDIM,NDIM)	
	DIMENSION ZBOGO1(NDIM,NDIM)
	DIMENSION ZBOGO2(NDIM,NDIM)
	
       
	
	call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZU1bar,NDIM,ZU0,NDIM,zzero,ZAUX1,NDIM)	
        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZV1bar,NDIM,ZV0,NDIM,zzero,ZAUX2,NDIM)

!        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,U0,NDIM,V0,NDIM,0.d0,AUX3,NDIM)
!        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,V0,NDIM,U0,NDIM,0.d0,AUX4,NDIM)
!        call DGEMM ('n','t',NDIM,NDIM,NDIM,-1.d0,U0,NDIM,V0,NDIM,0.d0,AUX6,NDIM)

        do ii=1,NDIM
	  do jj=1,NDIM
	      ZBOGO1(ii,jj)=ZAUX1(ii,jj)+ZAUX2(ii,jj)
	  end do
	end do

!      ........ ZBOGO1 is not Hermitian

       call clingd(NDIM,NDIM,NDIM,NDIM,ZBOGO1,ZX,zd,ifl)
       znorm = cdsqrt(zd)

!        write(2,*),' Onishi: Checking the relation: U^+U + V^+V = I'
        do ii=1,NDIM
        do jj=1,NDIM
!        if(abs(ZBOGO1(ii,jj)).gt. 0.1) write(2,'(2i5,2f20.4)') ii,jj,ZBOGO1(ii,jj)
        end do
        end do

        return

        end

       subroutine Oishi_old(ZU0,ZV0,ZU1bar,ZV1bar,znorm,N)           
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))

      Dimension ZV0(N,N),ZU0(N,N)
      Dimension ZU1bar(N,N),ZV1bar(N,N)
      Dimension ZAUX1(N,N),zx(N,N)
       print *,'Onishi formula'

!      call ZGEMM ('c','n',N,N,N,zone,ZU1bar,N,ZU0,N,zzero,ZAXU1,N)
!      call ZGEMM ('c','n',N,N,N,zone,ZV1bar,N,ZV0,N,zone,ZAXU1,N)
!      do i=1,N
!      do j=1,N
!        ZAUX1(i,j) = zzero
!      do k=1,N
!        ZAUX1(i,j) = ZAUX1(i,j) + conjg(ZU1bar(k,i))*ZU0(k,j) + conjg(ZV1bar(k,i))*ZV0(k,j) 
!      enddo
!      enddo
!      enddo

      call clingd(N,N,N,N,ZAXU1,zx,zd,ifl)
      znorm = cdsqrt(zd) 

      return
      end
