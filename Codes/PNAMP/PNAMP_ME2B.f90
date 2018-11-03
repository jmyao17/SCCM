        SUBROUTINE ZJ2_0(ZRO,ZKAPA10,ZKAPA01,zJ2_MV,NDIM)
        use vaphfb_par
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZKAPA10(NDIM,NDIM)
        DIMENSION ZKAPA01(NDIM,NDIM)
        dimension zAUX1(NDIM,NDIM)
        dimension zAUX2(NDIM,NDIM)
        dimension zAUX3(NDIM,NDIM)
        dimension zAUX4(NDIM,NDIM)
        dimension zAUX5(NDIM,NDIM)
        dimension zAUX6(NDIM,NDIM)


        !X
!        print *, '...x...'

!      ......... Rho*J^2_x

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &              cME1B%ZJx2ME,NDIM,zro,NDIM,zzero,zAUX1,NDIM)
!      ......... Rho*J_x
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &              cME1B%ZJxME,NDIM,zro,NDIM,zzero,zAUX2,NDIM)

!      ......... kappa10*J_x
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &              cME1B%ZJxME,NDIM,ZKAPA10,NDIM,zzero,zAUX3,NDIM)

!      ......... kappa01*J^T_x
        call ZGEMM ('t','t',NDIM,NDIM,NDIM,zone,  &
     &              cME1B%ZJxME,NDIM,ZKAPA01,NDIM,zzero,zAUX4,NDIM)

!      ......... 5: (J_x*Rho)*(J_x*Rho)
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &              zAUX2,NDIM,zAUX2,NDIM,zzero,zAUX5,NDIM)

!      ......... 6: (J_x*kappa10)*(J^T_x*kappa01)
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &              zAUX3,NDIM,zAUX4,NDIM,zzero,zAUX6,NDIM)


        ztrx1=zzero
        ztrx2=zzero
        ztrx3=zzero
        ztrx4=zzero
        do ii=1,NDIM
         ztrx1=ztrx1+zAUX1(ii,ii)   ! sum J^2_x*Rho
         ztrx2=ztrx2+zAUX2(ii,ii)   ! sum J_x*Rho
         ztrx3=ztrx3+zAUX5(ii,ii)   ! sum (J_x*Rho)*(J_x*Rho)
         ztrx4=ztrx4+zAUX6(ii,ii)   ! sum (J_x*kappa10)*(J^T_x*kappa01)
        enddo

        zxterm=ztrx1+ztrx2*ztrx2-ztrx3+ztrx4  ! <J^2_x>

!        print *, 'trx1=',ztrx1
!        print *, 'trx2=',ztrx2*ztrx2
!        print *, 'trx3=',-ztrx3
!        print *, 'trx4=',ztrx4
!        print *, '<J^2_x>=',zxterm

!        print *, '...y...'

        !y
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,   &
     &  cME1B%ZJy2ME,NDIM,zro,NDIM,zzero,zAUX1,NDIM)
!     &  cME1B%ZJy2ME,NDIM,zro,NDIM,zzero,zAUX1y,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,   &
     &  cME1B%ZJyME,NDIM,zro,NDIM,zzero,zAUX2,NDIM)
!     &  cME1B%ZJyME,NDIM,zro,NDIM,zzero,zAUX2y,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,   &
     &  cME1B%ZJyME,NDIM,ZKAPA10,NDIM,zzero,zAUX3,NDIM)
!     &  cME1B%ZJyME,NDIM,ZKAPA10,NDIM,zzero,zAUX3y,NDIM)

        call ZGEMM ('t','t',NDIM,NDIM,NDIM,zone,   &
     &  cME1B%ZJyME,NDIM,ZKAPA01,NDIM,zzero,zAUX4,NDIM)
!     &  cME1B%ZJyME,NDIM,ZKAPA01,NDIM,zzero,zAUX4y,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,   &
     &  zAUX2,NDIM,zAUX2,NDIM,zzero,zAUX5,NDIM)
!     &  zAUX2y,NDIM,zAUX2y,NDIM,zzero,zAUX5y,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,   &
     &  zAUX3,NDIM,zAUX4,NDIM,zzero,zAUX6,NDIM)
!     &  zAUX3y,NDIM,zAUX4y,NDIM,zzero,zAUX6y,NDIM)


        ztry1=zzero
        ztry2=zzero
        ztry3=zzero
        ztry4=zzero
        do ii=1,NDIM
         ztry1=ztry1+zAUX1(ii,ii)
         ztry2=ztry2+zAUX2(ii,ii)
         ztry3=ztry3+zAUX5(ii,ii)
         ztry4=ztry4+zAUX6(ii,ii)
        enddo

        zyterm=ztry1+ztry2*ztry2-ztry3+ztry4


!        print *, '...z...'
        !z
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &  cME1B%ZJz2ME,NDIM,zro,NDIM,zzero,zAUX1,NDIM)
!     &  cME1B%ZJz2ME,NDIM,zro,NDIM,zzero,zAUX1z,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &  cME1B%ZJzME,NDIM,zro,NDIM,zzero,zAUX2,NDIM)
!     &  cME1B%ZJzME,NDIM,zro,NDIM,zzero,zAUX2z,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &  cME1B%ZJzME,NDIM,ZKAPA10,NDIM,zzero,zAUX3,NDIM)
!     &  cME1B%ZJzME,NDIM,ZKAPA10,NDIM,zzero,zAUX3z,NDIM)

        call ZGEMM ('t','t',NDIM,NDIM,NDIM,zone,  &
     &  cME1B%ZJzME,NDIM,ZKAPA01,NDIM,zzero,zAUX4,NDIM)
!     &  cME1B%ZJzME,NDIM,ZKAPA01,NDIM,zzero,zAUX4z,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &  zAUX2,NDIM,zAUX2,NDIM,zzero,zAUX5,NDIM)
!     &  zAUX2z,NDIM,zAUX2z,NDIM,zzero,zAUX5z,NDIM)

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &  zAUX3,NDIM,zAUX4,NDIM,zzero,zAUX6,NDIM)
!     &  zAUX3z,NDIM,zAUX4z,NDIM,zzero,zAUX6z,NDIM)


        ztrz1=zzero
        ztrz2=zzero
        ztrz3=zzero
        ztrz4=zzero
        do ii=1,NDIM
         ztrz1=ztrz1+zAUX1(ii,ii)
         ztrz2=ztrz2+zAUX2(ii,ii)
         ztrz3=ztrz3+zAUX5(ii,ii)
         ztrz4=ztrz4+zAUX6(ii,ii)
        end do

        zzterm=ztrz1+ztrz2*ztrz2-ztrz3+ztrz4


        zJ2_MV=zxterm+zyterm+zzterm

!        print *, '<J^2_x>=', zxterm
!        print *, '<J^2_y>=', zyterm
!        print *, '<J^2_z>=', zzterm
        end subroutine
