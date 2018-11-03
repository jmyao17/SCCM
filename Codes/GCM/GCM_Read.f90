      subroutine Reader(lpr)
!    ...................
!    eta1     -> eta
!    Nitermax - > itermax
!    ...................

       USE VAPHFB_Par
       implicit none 
       logical lpr
       character*10 ctemp,a(50),b(50)
       integer neps,ii,isum
!       ................... read parameters       
        open(5, file='input.dat', status='old')
        read(5,'(A10,i1)')    b(1),Input%IntType
        read(5,'(A10,i1)')    ctemp,Input%IsHFB
!       .................... print out
        if(lpr) then
        write(*,*) ' ................................................. '
        print *,   ' start to read the input file'
        write(*,*) ' ................................................. '
        write(6,'(A)')   '....................'
        write(6,'(A)')   ' Input Parameters   '
        write(6,'(A)')   '....................'
        write(6,'(A10,i1)')   b(1),Input%IntType
        endif


!       ................. wfs by shell model NN int
        if(Input%IntType.eq.0) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20):: Input%cIntID)
           read(5,'(A10,a)')  b(2),Input%cIntID

           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6):: Input%cValID)
           read(5,'(A10,a6)')   b(3),Input%cValID
           read(5,'(A10,i2)')  ctemp,Input%ihwHO  ! newly added
           write(*,'(a,a20)') 'Interaction Type: ',Input%cIntID
           write(*,'(a,a6)') 'Model Space: ',trim(Input%cValID)
         endif
!       ................. wfs by chiral NN int
         if(Input%IntType.eq.1) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=8)  :: Input%cValID)
           read(5,'(A10,a20)')  b(2),Input%cIntID
           read(5,'(A10,a8)')   b(3),Input%cValID
           read(5,'(A10,i2)')  ctemp,Input%ihwHO  ! newly added
           if(lpr) write(6,'(A10,a20)')  b(2),Input%cIntID
           if(lpr) write(6,'(A10,a8)')   b(3),Input%cValID
           write(*,'(A10,a8)')   b(3),Input%cValID
         endif
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'GCN') &
     &  stop ' Interation Type is Wrong !'
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'KB3') &
     &  stop ' Interation Type is Wrong !'
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'USD') &
     &  stop ' Interation Type is Wrong !'

        if(.NOT. ALLOCATED(Input%cFlow)) ALLOCATE(character(len=4) :: Input%cFlow)
        read(5,'(A10,a4)')   b(5),Input%cFlow
!       ..................READ THE INPUT PARAMETERS
        read(5,'(A10,I6)')   b(7),Input%nprot
        read(5,'(A10,I6)')   b(8),Input%nneut
        read(5,'(A10,I6)')   b(9),PNP%NFOM
        PNP%MPhi = PNP%NFOM
        read(5,'(A10,I6)')   b(10),AMP%NLEG_ALP
        read(5,'(A10,I6)')   b(11),AMP%NLEG_BET
        read(5,'(A10,I6)')   b(12),AMP%NLEG_GAM
        read(5,'(A10,I6)')   b(13),GCM%Jmax
        if(GCM%Jmax.gt.JJmax) &
     &  stop 'Error: JJmax parameter should be larger !'
! .............................
        if(PNP%NFOM .le. 1)   PNP%NFOM = 1
        if(AMP%NLEG_ALP*AMP%NLEG_BET*AMP%NLEG_GAM .ge. 1)   AMP%IsAMP = 1
        if(AMP%NLEG_ALP*AMP%NLEG_GAM.gt.1) then
           AMP%i3DAMP = 1
        else
           AMP%i3DAMP = 0
        endif

        read(5,'(A10,I6)')    b(21),GCM%kmax
        if(GCM%kmax.gt.10)    stop ' GCM%kmax should be smaller than 10'
        read(5,'(A10,I6)')    b(18),Input%idens
        read(5,'(A10,I6)')    b(19),Input%iRho3B
        read(5,'(A10,a)')     ctemp,Input%vs4me3b
        if(Input%iRho3B .eq.1) print *, ' Model_Space for Rho3B:',Input%vs4me3b

!        read(5,'(A10,I6)')    b(22),Input%IE_Conv
        read(5,'(a10,1e9.2)')  b(20),GCM%eps
        print *, ' cutoff in norm:',GCM%eps
        Input%NOSJ(:) = 0     ! initialization
        read(5,'(A10,I6)')    b(22),Input%J0k   ! print out the w.f. of the 0^+_k  state
        read(5,'(A10,I6)')    b(23),Input%NOSJ(0)
        read(5,'(A10,I6)')    b(24),Input%NOSJ(2)
        read(5,'(A10,I6)')    b(24),Input%NOSJ(4)
        read(5,'(A10,I6)')    b(25),Input%icr
!  ............................

        ACONS_MV(1)= Input%nprot
        ACONS_MV(2)= Input%nneut
!   ................
        return
       end
