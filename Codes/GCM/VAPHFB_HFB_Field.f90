!    ...................................
!    From interaction to HFB fileds
!    ...................................
        subroutine HFB_FIELD_COMPLEX(zro,zkapa10,zkapa01,   &
     &  zakin,zgamma,zham,zdelta10,zdelta01,NLEV)
!      ..............................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)

                
        DIMENSION zro(NLEV,NLEV)
        DIMENSION zkapa10(NLEV,NLEV)
        DIMENSION zkapa01(NLEV,NLEV)
        DIMENSION zakin(NLEV,NLEV)
        DIMENSION zgamma(NLEV,NLEV)
        DIMENSION zham(NLEV,NLEV)
        DIMENSION zdelta10(NLEV,NLEV)
        DIMENSION zdelta01(NLEV,NLEV)
!        DIMENSION V_ABCD(NLEV,NLEV,NLEV,NLEV)

!     ...................... initialization
      do ia=1,NLEV
      do ic=1,NLEV
            ZGamma(ia,ic)  =zzero
            zDelta10(ia,ic)=zzero
            zDelta01(ia,ic)=zzero
            zham(ia,ic)    =zzero
      enddo
      enddo

      do iabcd=1,H%iabcd_max
          ia = H%ka(iabcd)
          ib = H%kb(iabcd)
          ic = H%kc(iabcd)
          id = H%kd(iabcd)
          ZGamma(ia,ic)=ZGamma(ia,ic) + H%ME2BM(iabcd)*ZRO(id,ib)  
          ZGamma(ib,ic)=ZGamma(ib,ic) - H%ME2BM(iabcd)*ZRO(id,ia)  
!          if(abs(H%ME2BM(iabcd)).gt.1.d-6)  write(*,*) ia,ib,ic,id,H%ME2BM(iabcd)
       enddo ! iabcd
!  ................... single-particle h 
       do ia=1,NLEV
         do ic=1,NLEV
            zham(ia,ic)=zGamma(ia,ic)+zakin(ia,ic)
      end do
        end do
        !Delta
      do iabcd=1,H%iabcd_max
          ia = H%ka(iabcd)
          ib = H%kb(iabcd)
          ic = H%kc(iabcd)
          id = H%kd(iabcd)
            zDelta10(ia,ib)=zDelta10(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa10(ic,id)
            if(ia.ne.ib) &
     &      zDelta10(ib,ia)=zDelta10(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa10(ic,id)

            zDelta01(ia,ib)=zDelta01(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa01(ic,id)
            if(ia.ne.ib) &
     &      zDelta01(ib,ia)=zDelta01(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa01(ic,id)

       enddo ! iabcd 
         
        end subroutine
