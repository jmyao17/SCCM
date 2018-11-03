!    ...................................
!       From interaction to HFB fileds
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
        DIMENSION demax(0:1)
!     ...................... initialization
      do ia=1,NLEV
      do ic=1,NLEV
            ZGamma(ia,ic)  =zzero
            zDelta10(ia,ic)=zzero
            zDelta01(ia,ic)=zzero
            zham(ia,ic)    =zzero
      enddo
      enddo

      HFB%ide(0) = 1      
      HFB%ide(1) = 1      
      do iabcd=1,H%iabcd_max

!       ................. introduce a cutoff in ME2B(M)
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

       enddo ! iabcd

!  ................... single-particle h 
       do ia=1,NLEV
         do ic=1,NLEV
            zham(ia,ic)=zGamma(ia,ic)+zakin(ia,ic) &
                        - shared%cf*cME1B%AJX_ME(ic+(ia-1)*NLEV)  ! cranking devided by 2
!  ..................w/o np-mixing and pairing
            if(Input%NPMix .eq. 0 .and. tnljm%t(ia) .ne. tnljm%t(ic)) then 
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif 
!  ..................w/o parity-mixing and pairing
            if(Input%IPMix .eq. 0 .and. iv(tnljm%l(ia)) .ne. iv(tnljm%l(ic))) then
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif
!  ..................w/o k-mixing (triaxiality) and pairing
            if(Input%KMix .eq. 0 .and. abs(tnljm%twom(ia)) .ne. abs(tnljm%twom(ic))) then
               zham(ia,ic)  = 0.d0
               zGamma(ia,ic)= 0.d0
            endif
         end do
        end do



!     HF/HFB
      if(Input%IsHFB.eq.0) then
       HFB%ide(0) = 0
       HFB%ide(1) = 0
       return
      endif
!      ............ Delta
      do iabcd=1,H%iabcd_max
          ia = H%ka(iabcd)
          ib = H%kb(iabcd)
          ic = H%kc(iabcd)
          id = H%kd(iabcd)
            zDelta10(ia,ib)=zDelta10(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa10(ic,id)*2
          
            if(ia.ne.ib) zDelta10(ib,ia)=zDelta10(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa10(ic,id)*2

            zDelta01(ia,ib)=zDelta01(ia,ib)+0.5*H%ME2BM(iabcd)*zkapa01(ic,id)*2

            if(ia.ne.ib) zDelta01(ib,ia)=zDelta01(ib,ia)-0.5*H%ME2BM(iabcd)*zkapa01(ic,id)*2

          if(ia+ib.eq.ic+id) cycle
!         .... (ab) <-> (cd)
            zDelta10(ic,id)=zDelta10(ic,id)+0.5*H%ME2BM(iabcd)*zkapa10(ia,ib)*2
            zDelta10(id,ic)=zDelta10(id,ic)-0.5*H%ME2BM(iabcd)*zkapa10(ia,ib)*2
            zDelta01(ic,id)=zDelta01(ic,id)+0.5*H%ME2BM(iabcd)*zkapa01(ia,ib)*2
            zDelta01(id,ic)=zDelta01(id,ic)-0.5*H%ME2BM(iabcd)*zkapa01(ia,ib)*2
       enddo ! iabcd

!  ..................w/o np-mixing
       do ia=1,NLEV
         do ib=1,NLEV
            if(Input%NPMix .eq. 0 .and. tnljm%t(ia) .ne. tnljm%t(ib)) then
             zDelta10(ia,ib) = 0.d0
             zDelta01(ia,ib) = 0.d0
            endif

            if(Input%IPMix .eq. 0 .and. iv(tnljm%l(ia)) .ne. iv(tnljm%l(ib))) then
             zDelta10(ia,ib) = 0.d0
             zDelta01(ia,ib) = 0.d0
            endif

            if(Input%KMix .eq. 0 .and. abs(tnljm%twom(ia)) .ne. abs(tnljm%twom(ib)) ) then
             zDelta10(ia,ib) = 0.d0
             zDelta01(ia,ib) = 0.d0
            endif

       enddo ! ib 
       enddo ! ia

!     ................. pairing collaps or not?
       demax(0:1) = zero
       do it=0,1
       do ia=it*NLEV/2+1,it*NLEV/2+NLEV/2
          do ib=it*NLEV/2+1,it*NLEV/2+NLEV/2
             if(abs(zDelta10(ia,ib)).gt. demax(it)) demax(it) = abs(zDelta10(ia,ib))
          enddo
       enddo
       if(abs(demax(it)).lt.0.005) then
           HFB%ide(it) = 0      ! label HFB wavefunction as HF wavefunction
           if(it.eq.0) print *, ' pairing collaps between neutrons ...'
           if(it.eq.1) print *, ' pairing collaps between protons ...'
       endif
       enddo


!   ........ pairing collaps for only either neutron or proton
       if(HFB%ide(it)+HFB%ide(it).eq.1) then
       do it=0,1
          if(HFB%ide(it).ne.0) cycle
          do ia=it*NLEV/2+1,it*NLEV/2+NLEV/2
          do ib=it*NLEV/2+1,it*NLEV/2+NLEV/2
             zDelta10(ia,ib) = 0.1 
             zDelta10(ib,ia) =-0.1 
             zDelta01(ia,ib) =-0.1 
             zDelta01(ib,ia) = 0.1 
          enddo
          enddo
       enddo
       endif


        end subroutine
