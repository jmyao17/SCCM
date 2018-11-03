           subroutine Model_Space()
           USE VAPHFB_PAR
           implicit none

           integer lj,nlj,nlj1,nlj2,t,n,l,twoj

           integer ie
           
           !write(*,*) ' Begin Model_Space: '
           ie=31
!          .......................... allocate

           if(.NOT. ALLOCATED(SPB%n)) ALLOCATE(SPB%n(1:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%l)) ALLOCATE(SPB%l(1:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%lj)) ALLOCATE(SPB%lj(1:HO%nljmax))
           if(.NOT. ALLOCATED(SPB%twoj)) ALLOCATE(SPB%twoj(1:HO%nljmax))
           HO%NOrbit(0) = 0
           HO%NOrbit(1) = 0
           open(ie, file="../val/"//trim(Input%cValID)//".val",&
                    status="old")
!          .......................... 

           read(ie,*)
           read(ie,*)
           read(ie,*)

           if(.NOT. ALLOCATED(SPB%nlj))&
          &ALLOCATE(SPB%nlj(0:HO%nmax,0:HO%ljmax))
           do nlj = 1, HO%nljmax
              read(ie,*) nlj1, t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
          if(nlj1.ne.nlj) stop ' Error in the definition of model space'
              lj = (SPB%l(nlj)*2 + SPB%twoj(nlj) -1)/2
              SPB%nlj(SPB%n(nlj),lj) = nlj

              SPB%lj(nlj) = lj
              HO%NOrbit(0) = HO%NOrbit(0) + (SPB%twoj(nlj)+1)
!              write(*,'(4i6)') t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
           enddo
!           write(*,'(a,i5)') '   #  of levels for neutron:',HO%NOrbit(0)

           read(ie,*)
!           write(*,*) '--------------------------------'
           do nlj = 1, HO%nljmax
              read(ie,*) nlj2, t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
          if(nlj2.ne.nlj) stop ' Error in the definition of model space'
              HO%NOrbit(1) = HO%NOrbit(1) + (SPB%twoj(nlj)+1)
!              write(*,'(4i6)') t, SPB%n(nlj),SPB%l(nlj),SPB%twoj(nlj)
           enddo
!           write(*,'(a,i5)') '   #  of levels for proton:',HO%NOrbit(1)
!           write(*,*) '--------------------------------'


          if(HO%NLEV .ne. HO%NOrbit(0)+HO%NOrbit(1)) then
            write(*,*) 'NLEV should be',HO%NOrbit(0)+HO%NOrbit(1)
            stop
          endif
          ! write(*,*) ' End Model_Space: '
           
         end subroutine Model_Space
