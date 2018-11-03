           subroutine Model_Space()
           USE VAPHFB_PAR
           implicit none

           integer lj,nlj,nlj1,nlj2,t,n,l,twoj

!          .......................... allocate

           if(.NOT. ALLOCATED(SPB%n)) ALLOCATE(SPB%n(0:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%l)) ALLOCATE(SPB%l(0:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%lj)) ALLOCATE(SPB%lj(0:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%twoj)) ALLOCATE(SPB%twoj(0:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%nlj))  ALLOCATE(SPB%nlj(0:HO%nmax,0:HO%ljmax))
!           write(*,*) '--------------------------------'
          do nlj=1,HO%nljmax
              SPB%n(nlj)    = -1
              SPB%l(nlj)    = -1
              SPB%lj(nlj)   = -1
              SPB%twoj(nlj) = -1
          enddo

           HO%NOrbit(0) = 0
           HO%NOrbit(1) = 0
           if(.NOT. ALLOCATED(SPB%VType)) &
      &    ALLOCATE(SPB%VType(0:HO%ljmax,0:HO%nmax,0:1))

           SPB%VType(0:HO%ljmax,0:HO%nmax,0:1) = -1

           open(21, file="../val/"//trim(Input%cValID)//".val",status="old")

!          .......................... 

           read(21,*)
           read(21,*)
           read(21,*)

           do nlj = 1, HO%nljmax
              read(21,*) nlj1, t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
              if(nlj1.ne.nlj) &
     &        stop ' Error in the definition of model space'
              lj = (SPB%l(nlj)*2 + SPB%twoj(nlj) -1)/2
              SPB%nlj(SPB%n(nlj),lj) = nlj
              SPB%lj(nlj) = lj
              HO%NOrbit(0) = HO%NOrbit(0) + (SPB%twoj(nlj)+1)
!              write(*,'(4i6)') t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
              SPB%VType(lj,SPB%n(nlj),t)   = 1
           enddo
!           write(*,'(a,i5)') '   #  of levels for neutron:',HO%NOrbit(0)

           read(21,*)
!           write(*,*) '--------------------------------'
           do nlj = 1, HO%nljmax
              read(21,*) nlj2, t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
              if(nlj2.ne.nlj) &
     &        stop ' Error in the definition of model space'
              HO%NOrbit(1) = HO%NOrbit(1) + (SPB%twoj(nlj)+1)
!              write(*,'(4i6)') t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
              lj = (SPB%l(nlj)*2 + SPB%twoj(nlj) -1)/2
              SPB%VType(lj,SPB%n(nlj),t)   = 1
           enddo
!           write(*,'(a,i5)') '   #  of levels for proton:',HO%NOrbit(1)
           write(*,*) '--------------------------------'


          if(HO%NLEV .ne. HO%NOrbit(0)+HO%NOrbit(1)) then
            write(*,*) 'NLEV should be',HO%NOrbit(0)+HO%NOrbit(1)
            stop
          endif

         end subroutine Model_Space
