!    ..............................................
!    From Hamiltonian (interaction) to HFB fileds
!    .............................................
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

          ZGamma(ia,id)=ZGamma(ia,id) - H%ME2BM(iabcd)*ZRO(ic,ib)
          ZGamma(ib,id)=ZGamma(ib,id) + H%ME2BM(iabcd)*ZRO(ic,ia)

          if(ia+ib.eq.ic+id) cycle

!         .... (ab) <-> (cd)
          ZGamma(ic,ia)=ZGamma(ic,ia) + H%ME2BM(iabcd)*ZRO(ib,id)
          ZGamma(id,ia)=ZGamma(id,ia) - H%ME2BM(iabcd)*ZRO(ib,ic)

          ZGamma(ic,ib)=ZGamma(ic,ib) - H%ME2BM(iabcd)*ZRO(ia,id)
          ZGamma(id,ib)=ZGamma(id,ib) + H%ME2BM(iabcd)*ZRO(ia,ic)


!          if(abs(H%ME2BM(iabcd)).gt.1.d-6) write(201,*) ia,ib,ic,id,H%ME2BM(iabcd)
       enddo ! iabcd

!  ................... single-particle h 
       do ia=1,NLEV
         do ic=1,NLEV
            zham(ia,ic)=zGamma(ia,ic)+zakin(ia,ic)
!  ..................w/o np-mixing and pairing
            if(Input%NPMix .eq. 0 .and. tnljm%t(ia) .ne. tnljm%t(ic)) then 
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif 
!  ..................w/o np-mixing and pairing
            if(Input%IPMix .eq. 0 .and. iv(tnljm%l(ia)) .ne. iv(tnljm%l(ic))) then
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif
!  ..................w/o k-mixing (triaxiality) and pairing
            if(Input%KMix .eq. 0 .and. tnljm%twom(ia) .ne. tnljm%twom(ic)) then
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif
         end do
        end do
        !Delta
      do iabcd=1,H%iabcd_max
          ia = H%ka(iabcd)
          ib = H%kb(iabcd)
          ic = H%kc(iabcd)
          id = H%kd(iabcd)
            zDelta10(ia,ib)=zDelta10(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa10(ic,id)*2

            zDelta10(ib,ia)=zDelta10(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa10(ic,id)*2

            zDelta01(ia,ib)=zDelta01(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa01(ic,id)*2

            zDelta01(ib,ia)=zDelta01(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa01(ic,id)*2

          if(ia+ib.eq.ic+id) cycle
!         .... (ab) <-> (cd)
            zDelta10(ic,id)=zDelta10(ic,id)+0.5*H%ME2BM(iabcd)*zkapa10(ia,ib)*2
            zDelta10(id,ic)=zDelta10(id,ic)-0.5*H%ME2BM(iabcd)*zkapa10(ia,ib)*2
            zDelta01(ic,id)=zDelta01(ic,id)+0.5*H%ME2BM(iabcd)*zkapa01(ia,ib)*2
            zDelta01(id,ic)=zDelta01(id,ic)-0.5*H%ME2BM(iabcd)*zkapa01(ia,ib)*2
       enddo ! iabcd 
!  ..................w/o np-mixing and pairing
       do ia=1,NLEV
         do ib=1,NLEV
            if(Input%NPMix .eq. 0 .and. tnljm%t(ia) .ne. tnljm%t(ib)) then 
             zDelta10(ia,ib) = 0.d0
             zDelta10(ib,ia) = 0.d0
             zDelta01(ia,ib) = 0.d0
             zDelta01(ib,ia) = 0.d0
            endif
            if(Input%IPMix .eq. 0 .and. iv(tnljm%l(ia)) .ne. iv(tnljm%l(ib))) then
             zDelta10(ia,ib) = 0.d0
             zDelta10(ib,ia) = 0.d0
             zDelta01(ia,ib) = 0.d0
             zDelta01(ib,ia) = 0.d0
            endif

            if(Input%KMix .eq. 0 .and. abs(tnljm%twom(ia)) .ne. abs(tnljm%twom(ib)) ) then
             zDelta10(ia,ib) = 0.d0
             zDelta10(ib,ia) = 0.d0
             zDelta01(ia,ib) = 0.d0
             zDelta01(ib,ia) = 0.d0
            endif

       enddo ! ib 
       enddo ! ia 
         
        end subroutine
