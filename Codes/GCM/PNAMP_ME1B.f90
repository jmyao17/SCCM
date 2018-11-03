!      .................................................
        SUBROUTINE ZQ2(ZQME,ZRO,zQ_0p,zQ_0n,NDIM)
!      .................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)


        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)
        zQ_0p=zzero
        zQ_0n=zzero
        do ii=1,NDIM/2
         zQ_0p=zQ_0p+ZAUX(ii,ii)
         zQ_0n=zQ_0n+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE

!      .................................................
        SUBROUTINE ZQ2t(ZQME,ZRO,zQ_0p,zQ_0n,NDIM)
!      .................................................
!       transponent of Q2
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)


        call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)
        zQ_0p=zzero
        zQ_0n=zzero
        do ii=1,NDIM/2
         zQ_0p=zQ_0p+ZAUX(ii,ii)
         zQ_0n=zQ_0n+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE
