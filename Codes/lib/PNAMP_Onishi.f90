      subroutine Onishi(ZU0,ZV0,ZU1bar,ZV1bar,znorm,NDIM)           
!    ..................
!    Onishi formula
!    ..................

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

!     attention: not 'c', but 't'       

        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZU1bar,NDIM,ZU0,NDIM,zzero,ZAUX1,NDIM)	
        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,ZV1bar,NDIM,ZV0,NDIM,zzero,ZAUX2,NDIM)

        do ii=1,NDIM
          do jj=1,NDIM
             ZBOGO1(ii,jj)=ZAUX1(ii,jj)+ZAUX2(ii,jj)
          end do
        end do

!      ........ ZBOGO1 is not Hermitian
       ! calculate the determinant of the U
       call clingd(NDIM,NDIM,NDIM,NDIM,ZBOGO1,ZX,zd,ifl)
       znorm = cdsqrt(zd)

!        write(2,*),' Onishi: Checking the relation: U^+U + V^+V = I'
!        do ii=1,NDIM
!        do jj=1,NDIM
!        if(abs(ZBOGO1(ii,jj)).gt. 0.1) write(2,'(2i5,2f20.4)') ii,jj,ZBOGO1(ii,jj)
!        end do
!        end do

        return

        end

