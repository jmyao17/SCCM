         subroutine canon(NDIM)
         use VAPHFB_Par
         implicit none
         integer m1,m2,mm,ii,jj,NDIM,INFO2
         real*8  ez(NDIM),Ocu_rho(NDIM),roaux(NDIM,NDIM)
         real*8  WORKrho(3*NDIM-1)
         real*8  AUX11(NDIM,NDIM),AUX12(NDIM,NDIM)


       


       if(.NOT. allocated(HFB%it1)) allocate(HFB%it1(1:NDIM))
       if(.NOT. allocated(HFB%it2)) allocate(HFB%it2(1:NDIM))
       if(.NOT. allocated(HFB%ip1)) allocate(HFB%ip1(1:NDIM))
       if(.NOT. allocated(HFB%ip2)) allocate(HFB%ip2(1:NDIM))

       if(.NOT. allocated(HFB%RO_1)) &
     & allocate(HFB%RO_1(1:NDIM,1:NDIM))
       if(.NOT. allocated(HFB%vv1)) allocate(HFB%vv1(1:NDIM))

       if(.NOT. allocated(HFB%RO_2)) &
     & allocate(HFB%RO_2(1:NDIM,1:NDIM))
       if(.NOT. allocated(HFB%vv2)) allocate(HFB%vv2(1:NDIM))

!    ............. density matrix of state 1 (U_0,V_0) 
!    RO = V^* V^T
       call DGEMM ('n','t',NDIM,NDIM,NDIM,-1.d0,               &
     &             HFB%V0,NDIM,HFB%V0,NDIM,0.d0,HFB%RO_1,NDIM)
!  ................w/o np-mixing and pairing
      if(Input%NPMix .eq. 0) then
      do m1=1,NDIM
      do m2=1,NDIM
         if(tnljm%t(m1) .ne. tnljm%t(m2)) then
            HFB%RO_1(m1,m2)    = 0.d0
          endif
       enddo
       enddo
       endif

!       call test()
       call dsyev('V','U',NDIM,HFB%RO_1,NDIM,Ocu_rho,WORKrho,3*NDIM-1,INFO2)
!       call sdiag(ndim,ndim,HFB%RO_1,Ocu_rho,HFB%RO_1,ez,1)

!    Using  U_rho to transform t3 
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  HFB%RO_1,NDIM,t3me,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,HFB%RO_1,NDIM,zero,AUX12,NDIM)

!    ........... the neutron and proton states are mixed in the index jj
!    mm index: p1,p2,...,pNLEV/2, n1,n2,...,nNLEV/2
!    jj index: ordered in energy
!   ........................ 
!      do jj=1,NDIM
!      do mm=1,NDIM
!       if(abs(HFB%RO_1(mm,jj)).gt.1.d-3) write(*,'(2i3,1f8.5)') jj,mm, HFB%RO_1(mm,jj)
!       enddo
!      enddo

      do ii=1,NDIM
          HFB%it1(ii) =AUX12(ii,ii)
         if(-Ocu_rho(ii).ge.1.d0) then
           HFB%vv1(ii) =one
         elseif(-Ocu_rho(ii).le.0.d0) then
           HFB%vv1(ii) = zero
         else
          HFB%vv1(ii) =-Ocu_rho(ii)
         endif
      enddo 

!    ............. parity
!    Using  U_rho to transform ip 
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  HFB%RO_1,NDIM,ipme,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,HFB%RO_1,NDIM,zero,AUX12,NDIM) 

      do ii=1,NDIM
          HFB%ip1(ii) =AUX12(ii,ii)
      enddo 



!    ............. density matrix of state 2 (U_1,V_1)
!    RO = V^* V^T
       call DGEMM ('n','t',NDIM,NDIM,NDIM,-1.d0,               &
     &             HFB%V1,NDIM,HFB%V1,NDIM,0.d0,HFB%RO_2,NDIM)
!  ................w/o np-mixing and pairing
      if(Input%NPMix .eq. 0) then
      do ii=1,NDIM
      do jj=1,NDIM
         if(tnljm%t(ii) .ne. tnljm%t(jj)) then
            HFB%RO_2(ii,jj)    = 0.d0
          endif
       enddo
       enddo
       endif

       call dsyev('V','U',NDIM,HFB%RO_2,NDIM,Ocu_rho,WORKrho,3*NDIM-1,INFO2)
!       call sdiag(ndim,ndim,HFB%RO_2,Ocu_rho,HFB%RO_2,ez,1)
!    Using  U_rho to transform t3 

!      AUX12: t3(jj,kk) = sum_(m1,m2) RO_2(m1,jj) * t3(m1,m2)* RO_2(m2,kk)      

        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  HFB%RO_2,NDIM,t3me,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,HFB%RO_2,NDIM,zero,AUX12,NDIM)

!      do jj=1,NDIM
!      do mm=1,NDIM
!       if(abs(HFB%RO_2(mm,jj)).gt.1.d-3) write(*,'(2i3,1f8.5)') jj,mm, HFB%RO_2(mm,jj)
!       enddo
!      enddo

      do ii=1,NDIM
          HFB%it2(ii) =AUX12(ii,ii)
         if(-Ocu_rho(ii).ge.1.d0) then
           HFB%vv2(ii) = one
         elseif(-Ocu_rho(ii).le.0.d0) then
           HFB%vv2(ii) = zero
         else
          HFB%vv2(ii) =-Ocu_rho(ii)
         endif
!         write(*,*) ii,HFB%vv2(ii)
      enddo

!    ............. parity
!    Using  U_rho to transform ip 
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  HFB%RO_2,NDIM,ipme,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,HFB%RO_2,NDIM,zero,AUX12,NDIM) 
      
      do ii=1,NDIM
          HFB%ip2(ii) =AUX12(ii,ii)
      enddo
!   ............. print out    
!   ii: index for s.p. state, not ordered in isospin projection any more
!   ii is ordered in occupation probability.
      do ii=1,NDIM
         write(*,'(i3,4f12.7)') ii,HFB%vv1(ii),HFB%it1(ii),HFB%ip1(ii) !,HFB%vv2(ii),HFB%it2(ii)
      enddo
      stop
      end


!______________________________________________________________________________
       subroutine dimens(NDIM)
!     .........................................................................
!     determine the dimension of rotional matrix
!     .........................................................................
      use vaphfb_par
      implicit none
!      implicit real*8 (a-h,o-z)
      integer ii, npz1,npz2
      real*8 voc1(NDIM),voc2(NDIM)
      integer k,l,k1,k2,NDIM,l1,l2,ldi,kdi
      integer kdim,ldim

      if(.NOT. allocated(HFB%kocc1)) allocate(HFB%kocc1(1:NDIM))
      if(.NOT. allocated(HFB%kocc2)) allocate(HFB%kocc2(1:NDIM))

            k1 = 1       ! ka1(ib,it) + 1
            k2 = NDIM ! ka1(ib,it) + kd1(ib,it)
            kdi= NDIM ! kd1(ib,it)
            l1 = 1       !ka2(ib,it) + 1
            l2 = NDIM !ka2(ib,it) + kd2(ib,it)
            ldi= NDIM !kd2(ib,it)
!----- determine the dimension of rotational matrix
            do k=1,kdi
               voc1(k)=zero
            enddo ! k  
            do l=1,ldi
               voc2(l)=zero
            enddo ! k  
!----- determine the dimension of rotational matrix
           do k=k1,k2
              voc1(k-k1+1) = HFB%vv1(k)
           enddo !k
           do l=l1,l2
              voc2(l-l1+1) = HFB%vv2(l)
           enddo !l
!------ order the occupation probabilities
!           call ordls(kdi,voc1)
!           call ordls(ldi,voc2)

           kdim = 0
           do k=k1,k2
              if(HFB%vv1(k).le.PNP%eps) goto 10
              kdim = kdim +1
   10      enddo !k
              npz1=kdim
              kdim = 0
           do l=l1,l2
              if(HFB%vv2(l).le.PNP%eps) goto 20
              kdim = kdim +1
   20      enddo !l
              npz2=kdim
          if(npz1.lt.npz2) then
              PNP%npz     = npz1
              PNP%eps1= PNP%eps
              npz2    = npz1
              PNP%eps2= voc2(npz2+1)
           else if(npz1.gt.npz2) then
              PNP%npz      = npz2
              PNP%eps2 = PNP%eps
              npz1     = npz2
              PNP%eps1 = voc1(npz1+1)
           else
               PNP%eps1 = PNP%eps
               PNP%eps2 = PNP%eps
               PNP%npz      = npz1   ! or npz2(ib,it)
           endif
!    ...................................................... determine dimension of rotation matrix
           kdim = 0
           do k=k1,k2
              if(HFB%vv1(k).le.PNP%eps1) goto 50
                 kdim = kdim +1
                 HFB%kocc1(k)=kdim
   50      enddo !k  
              npz1=kdim
              ldim = 0
           do l=l1,l2
              if(HFB%vv2(l).le.PNP%eps2) goto 60
                 ldim = ldim +1
                 HFB%kocc2(l)=ldim
   60      enddo !l
              npz2=ldim
              if(npz1.ne.npz2) stop 'in dimens: npz1 neq npz2'
              PNP%npz = npz1

!      .....................
        write(*,*) 'npz=',PNP%npz
        write(*,*) 'eps1=',PNP%eps1,' eps2=',PNP%eps2
        do ii=1,PNP%npz
           write(*,*) ii,HFB%it1(ii),HFB%vv1(ii) !,HFB%HFB%vv2(ii)
       enddo ! ii
!      ...................
      return
!-END-dimension
      END

