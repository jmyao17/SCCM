!_______________________________________________________________________
       program HWG_GCM
!-----------------------------------------------------------------------  
      use VAPHFB_Par
      implicit none
      integer nmaxdj,jproj,iq1,iq2
      integer kmax,k1,k2,k1m,k2m,iki0
      integer k,iqk,iqk1,iqk2
      integer maxmp
      real*8, dimension(:,:), allocatable ::  FF
      real*8, dimension(:,:,:), allocatable :: Rho2BJ0
      real*8, dimension(:,:,:,:), allocatable :: Rho1BJ0
      real*8, DIMENSION(:,:,:,:,:,:), allocatable :: OBD_Save
      real*8, DIMENSION(:,:), allocatable :: Rho3BJ0  
      real*8 wignei
!.........................................................................
   99 format (/,'  __________________________________________________ ', &
     &        /,' |                                                  |', &
     &        /,' |  program gcm                      version 01     |', &
     &        /,' |                                                  |', &
     &        /,' |  01 Nov. 2017                                    |', &
     &        /,' |  Copyright @ J. M. Yao                           |', &
     &        /,' |                                                  |', &
     &        /,' |  If you have found any problem in the code       |', &
     &        /,' |                                                  |', &
     &        /,' |  please contact us by sending email.             |', &
     &        /,' |__________________________________________________|', &
     &        /)

  100 format ('      No Projections ')
  101 format(a,i2)
  102 format(/,1x,74('-'),/)
  103 format ('      PNP  ')
  104 format ('      AMP  ')
  200 format ('      PNP+AMP  ')
  201 format ('  J  K1 K2 iq     beta   gamma     qtot       P00       E &
      &      ',           '   n^J(q,q)    <N>     <Z>     <J^2>')

  400 format(/,' existing matrix elements',                         &
     &         ' -- value 99 means missing, otherwise ln(overlap)',/)
  401 format(2x,2i3,1x,35i3)
  402 format ('    K  iq')
  403 format ('                                  ')
      write(*,99)

!     ........................................................ default value
       call PNAMP_Default()
!     ........................................................ read parameters
       call Reader(.false.)
       call Prep(.false.)
       call ReadCoordinators()
       GCM%NOQ=Input%NGCM 
       call Model_space(.false.)
       call ZCONST
       call gfv
       call ME2B_Base_Full(.false.)
       if(Input%iRho3B.eq.1) then
        call Rho3B_1B_Model_Space()
        call Rho3B_2B_Model_Space()
        call ME3B_Base_Full(.true.)
       endif

      if(AMP%NLEG_Bet.eq.1.and.PNP%NFOM.eq.1) write(*,100)
      if(AMP%NLEG_Bet.eq.1.and.PNP%NFOM.gt.1) write(*,103)
      if(AMP%NLEG_Bet.gt.1.and.PNP%NFOM.eq.1) write(*,104)
      if(AMP%NLEG_Bet.gt.1.and.PNP%NFOM.gt.1) write(*,200)

      GCM%Jmin  = 0
      if(AMP%i3DAMP.eq.0) then
        GCM%Jdf  = 2
      elseif(AMP%i3DAMP.eq.1) then
        GCM%Jdf  = 1
      else
        stop ' AMP%i3DAMP should be either 0 or 1'
      endif

!     .................................................. initialization
      do iq1=1,NOQmax*JJmax*NOQmax*JJmax !(GCM%NOQ*GCM%Jmax)**2 !jmax*maxq*jmax
         GCM%iexst(iq1,1) = 99
      enddo
!............................... read kernels 
      print *, 'Number of configurations:',GCM%NOQ
      call Filename4collwf(GCM%NOQ)
      do iq1=1, GCM%NOQ     
      do iq2=iq1, GCM%NOQ     
         call Filename4Kernels(AMP%i3DAMP,iq1,iq2)
             open(lou,file=File%elem,err=500)
             ! calculate BE2
             ! print *,File%TD1B
             ! print *,File%elem
             ! print *,AMP%i3DAMP
             if(Input%icr.eq.0) then
               ! no cranking, GCM%nmaxdi(iis)  = (iis/2+1)*GCM%NOQ or ...
               if(AMP%i3DAMP.eq.0) call readkernel_axial(iq1,iq2,GCM%Jmin,GCM%Jmax,GCM%Jdf)
               if(AMP%i3DAMP.ne.0) call readkernel_general(iq1,iq2,GCM%Jmin,GCM%Jmax,GCM%Jdf)
               call HN_Matrix(GCM%Jmin,GCM%Jmax,GCM%Jdf,iq1,iq2,AMP%i3DAMP,.false.)
             else  
               ! with cranking, GCM%nmaxdi(iis)  = GCM%NOQ*(2iis+1)
               call readkernel_general(iq1,iq2,GCM%Jmin,GCM%Jmax,GCM%Jdf)
               call HN_Matrix_general(GCM%Jmin,GCM%Jmax,GCM%Jdf,iq1,iq2,.false.)
             endif
 500         close(lou)

      enddo
      enddo
!    ... a bug was here
!      if(.not. allocated(FJkq)) allocate(FJkq(1:GCM%NOQ*JJmax,0:JJmax,GCM%kmax))
!    ......................................................... loop over angular momentum J
      if(Input%icr.eq.0) then
         call solution_triaxial()
         call spectrum(GCM%Jmin,GCM%Jmax,GCM%Jdf)
         call trans_by_me(GCM%Jmin,GCM%Jmax,GCM%Jdf)
      else 
         call solution_general()
         call spectrum_general(GCM%Jmin,GCM%Jmax,GCM%Jdf)
         call trans_by_me_general(GCM%Jmin,GCM%Jmax,GCM%Jdf)
      endif


      if(AMP%i3DAMP.ne.0) stop ' ..... END ....'
! up here: dimension should be changed starting from here
!     ............................................................ calculate observables
!    ..............................................
!     ............................................................ print spectrum  
!     call spect(GCM%Jmin,GCM%Jmax,GCM%Jdf)
!     .................................................... print table with trans. strengths

       !print *, ' --> start to calculate densities'

       if( .NOT. allocated(GCMDens%TD1BSum)) &
     & allocate(GCMDens%TD1BSum(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax, &
                                0:JJmax,0:JJmax,1:10,1:10,0:2))
       GCMDens%TD1BSum = zero

!     ............ one-body density
       if( .NOT. allocated(GCMDens%Rho1BSum)) &
     & allocate(GCMDens%Rho1BSum(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax))  ! it, lj, n1,n2
       GCMDens%Rho1BSum = zero     

       if( .NOT. allocated(Rho1BJ0)) &
     & allocate(Rho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax))           ! it, lj, n1,n2

!     ..................... the 2B, 3B density matrices for ground state
!       if(Input%idens.ne.1) &
!     & stop 'Warning: Densities are not calculated !' 

       if( .NOT. allocated(OBD_Save)) &
     & allocate(OBD_Save(0:1,1:HO%nMax,0:HO%LMax*2,0:1,1:HO%nMax,0:HO%LMax*2))
       !print *, ' one-body density is allocated ...' 
!     .............................
       if(Input%idens.eq.1) then 
         if( .NOT. allocated(GCMDens%Rho2BSum)) &
     &   allocate(GCMDens%Rho2BSum(-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax))

         if( .NOT. allocated(Rho2BJ0)) &
     &   allocate(Rho2BJ0(-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax))

         if(Input%iRho3B.eq.1) then
            if( .NOT. allocated(GCMDens%Rho3BSum)) &
     &      allocate(GCMDens%Rho3BSum(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))  
            if( .NOT. allocated(GCMDens%L3BSum)) &
     &      allocate(GCMDens%L3BSum(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))  
            if( .NOT. allocated(Rho3BJ0)) &
     &      allocate(Rho3BJ0(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))  
          endif
      endif
!    ................ initialization

      write(*,102)
      if(Input%idens.eq.1) then 
!    ........ initialization
      GCMDens%Rho2BSum = zero     
      if(Input%iRho3B.eq.1) GCMDens%Rho3BSum = zero     
      endif
!    ......... For J=0, GCM%nmaxdi(ji=0) = GCM%NOQ 
      do iq1=1, GCM%NOQ
      do iq2=iq1, GCM%NOQ
         call Filename4Kernels(AMP%i3DAMP,iq1,iq2)
!        ............. 1B
         !print *, File%Rho1B
         open(lou,file=File%Rho1B,err=700)
700      call Read_Rho1B_J0(lou,Rho1BJ0,.false.)
         call Add_Rho1B_J0(iq1,iq2,Rho1BJ0,.false.)
         close(lou)

         !print *, File%TD1B
         call Calc_TD1B(iq1,iq2,.false.)

         if(Input%idens.ne.1) cycle
!        ............. 2B
!         print *, File%Rho2B
         open(lou,file=File%Rho2B,err=800)
800      call Read_Rho2B_J0(lou,Rho2BJ0,.false.)
         call Add_Rho2B_J0(iq1,iq2,Rho2BJ0,.false.)
         close(lou)
!        ............. 3B
         if(Input%iRho3B.eq.1) then
         open(lou,file=File%Rho3B,err=900)
900      call Read_Rho3B_J0(lou,Rho3BJ0,.false.)
         call Add_Rho3B_J0(iq1,iq2,Rho3BJ0,.false.)
         close(lou)
         endif
       enddo
       enddo

       write(*,102)
       call Check_Rho1B(.false.)
       call Write_Rho1B4IMSRG(OBD_Save)
       call diag_Rho1B(.false.)
       !call Write_TD1B(OBD_Save)

       if(Input%idens.eq.1) then
           call Check_Rho2B(.false.)
           call Write_Rho2B(222,OBD_Save,GCMDens%Rho2BSum)
 
          if(Input%iRho3B.eq.1) then
            call Rho3B2Lambda3B(400,GCMDens%Rho1BSum, &
                GCMDens%Rho2BSum,GCMDens%Rho3BSum,GCMDens%L3BSum)
             call Check_Rho3B(.true.)
             call Write_L3B(lou,GCMDens%L3BSum,.true.)
          endif
       endif

       ! computing observables using one-body transition density
       write(*,102)
       if(.false.) then
          print *, ' --> re-calculate transitions', &
                         ' with one-body transition densities'
          call trans_by_td(GCM%Jmin,GCM%Jmax,GCM%Jdf)
       endif

       call Write_TD1B(.true.)

      end

