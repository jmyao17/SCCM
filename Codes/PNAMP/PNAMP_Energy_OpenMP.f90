!    ...................................
!    From interaction to HFB fileds
!    eMax08: 11-13 seconds
!    ...................................
        subroutine HFB_Energy_COMPLEX(zro,zkapa10,zkapa01,zakin,zEkin,zEHFB,NLEV,lpr)
!      ..............................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)
        logical lpr
        integer :: myID,OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        real*8    Delta_max                
        DIMENSION zro(NLEV,NLEV)
        DIMENSION zkapa10(NLEV,NLEV)
        DIMENSION zkapa01(NLEV,NLEV)
        DIMENSION zakin(NLEV,NLEV)
        DIMENSION zgamma(NLEV,NLEV)
        DIMENSION zham(NLEV,NLEV)
        DIMENSION zdelta10(NLEV,NLEV)
        DIMENSION zdelta01(NLEV,NLEV)

        real*8 omp_get_wtime,t0,t1
      if(lpr) then
        t0 = omp_get_wtime()
        write(*,'(a30)',advance='no') '  Computing the <HR>/<R> ...'
      endif
      zEHFB = zzero 
! Using the OpenMP takes even longer time
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(zEHFB,H,ZRO,ZKAPA10,ZKAPA01) 
!!$OMP DO SCHEDULE(DYNAMIC) REDUCTION(+:zEHFB)
      do iabcd=1,H%iabcd_max !iabcd_max
          k_i = H%ka(iabcd)
          k_j = H%kb(iabcd)
          k_k = H%kc(iabcd)
          k_l = H%kd(iabcd)

          V2B = H%ME2BM(iabcd)


       zrho_ijkl = (ZRO(k_k,k_i)*ZRO(k_l,k_j)                &
     &              -ZRO(k_l,k_i)*ZRO(k_k,k_j)                &
     &              +ZKAPA01(k_i,k_j)*ZKAPA10(k_k,k_l))

       zrho_jikl = (ZRO(k_k,k_j)*ZRO(k_l,k_i)                &
     &              -ZRO(k_l,k_j)*ZRO(k_k,k_i)                &
     &              +ZKAPA01(k_j,k_i)*ZKAPA10(k_k,k_l))

       zrho_jilk = (ZRO(k_l,k_j)*ZRO(k_k,k_i)                &
     &              -ZRO(k_k,k_j)*ZRO(k_l,k_i)                &
     &              +ZKAPA01(k_j,k_i)*ZKAPA10(k_l,k_k))

        zrho_ijlk = (ZRO(k_l,k_i)*ZRO(k_k,k_j)                &
     &              -ZRO(k_k,k_i)*ZRO(k_l,k_j)                &
     &              +ZKAPA01(k_i,k_j)*ZKAPA10(k_l,k_k))

!   -------- ij <--> kl
        zrho_klij = (ZRO(k_i,k_k)*ZRO(k_j,k_l)                &
     &              -ZRO(k_j,k_k)*ZRO(k_i,k_l)                &
     &              +ZKAPA01(k_k,k_l)*ZKAPA10(k_i,k_j))

        zrho_klji = (ZRO(k_j,k_k)*ZRO(k_i,k_l)                &
     &              -ZRO(k_i,k_k)*ZRO(k_j,k_l)                &
     &              +ZKAPA01(k_k,k_l)*ZKAPA10(k_j,k_i))

        zrho_lkij = (ZRO(k_i,k_l)*ZRO(k_j,k_k)                &
     &              -ZRO(k_j,k_l)*ZRO(k_i,k_k)                &
     &              +ZKAPA01(k_l,k_k)*ZKAPA10(k_i,k_j))

        zrho_lkji = (ZRO(k_j,k_l)*ZRO(k_i,k_k)                &
     &              -ZRO(k_i,k_l)*ZRO(k_j,k_k)                &
     &              +ZKAPA01(k_l,k_k)*ZKAPA10(k_j,k_i))

       d_ij_kl = 1.d0

       if((k_k+k_l).eq.(k_i+k_j)) d_ij_kl= 0.d0 

        zEHFB = zEHFB + 0.25*V2B*(zrho_ijkl+d_ij_kl*zrho_klij-zrho_jikl-d_ij_kl*zrho_klji   &
     &                           -zrho_ijlk-d_ij_kl*zrho_lkij+zrho_jilk+d_ij_kl*zrho_lkji) 

       enddo ! iabcd 
!!$OMP  END DO
!!$OMP  END PARALLEL

        zEkin = zzero
        do ia=1,NLEV
         do ic=1,NLEV
            zEkin = zEkin +  ZRO(ic,ia) * zakin(ia,ic) ! H%ME1BM(ia,ic) * ZRO(ia,ic)  !zakin(ia,ic)            
         end do
        end do
            zEHFB = zEHFB + zEkin            

        if(lpr) then
          write(*,'(a5)',advance='no') 'done'
          t1 = omp_get_wtime()
          write(*,'(a20,f10.5)',advance='no'), '...Time elapsed:',t1-t0
          write(*,'(a14,2f15.8)') ' <HR>/<R>=',zEHFB 
        endif
        return

        end subroutine
