      subroutine Reader(lpr)
!    ...................
!    eta1     -> eta
!    Nitermax - > itermax
!    ...................

       USE VAPHFB_Par
       implicit none 
       logical lpr
       character*10 ctemp,a(0:50),b(0:50)
       integer itemp,neps,ii,isum
       
       if(lpr) then
        print *, ' ======================== '
        print *, ' Features of Calculations '
        print *, ' ======================== '
       endif

! --- Open input files
        open(lin, file='input.dat', status='old')
        read(lin,'(A10,i1)')    b(0),Input%NPMix
        read(lin,'(A10,i1)')    ctemp,Input%IPMix
        read(lin,'(A10,i1)')    ctemp,Input%KMix
        read(lin,'(A10,i1)')    ctemp,Input%IsHFB
        read(lin,'(A10,i1)')    b(1),Input%IntType

       if(lpr) then
        if(Input%NPMix.eq.0) then
           print *, ' -> No n-p mixing! '
        else
           print *, ' -> Allow for n-p mixing ^_^ '
        endif
        if(Input%IPMix.eq.0) then
           print *, ' -> No parity mixing! '
        else
           print *, ' -> Allow for parity mixing ^_^ '
        endif
        if(Input%KMix.eq.0) then
           print *, ' -> No K mixing/Triaxiality! '
        else
           print *, ' -> Allow for K-mixing/Triaxiality ^_^ '
        endif

        if(Input%IntType.eq.0) then
           print *,' -> Shell-Model Interaction in Limited Model Space!'
        else
           print *,' -> Chiral Interaction in Full Model Space ^_^ '
        endif
        endif
!       ................. wfs by shell model NN int
!       ................. shell model int
        if(Input%IntType.eq.0) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           read(lin,'(A10,a20)')  b(2),Input%cIntID

           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6) :: Input%cValID)
           read(lin,'(A10,a5)')   b(3),Input%cValID

           read(lin,'(A10,i2)')  ctemp,Input%ihwHO
           read(lin,'(A10,i2)')  ctemp,Input%iCOM

           write(Input%chwHO,'(a4,i2)') 'hwHO',Input%ihwHO
!           write(*,'(a,a5)') '     Model Space: ',trim(Input%cValID)
         endif

!       ................. wfs by chiral NN int
         if(Input%IntType.eq.1) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6)  :: Input%cValID)
           read(lin,'(A10,a20)')  b(2),Input%cIntID
           read(lin,'(A10,a6)')   b(3),Input%cValID
           read(lin,'(A10,i2)')  ctemp,Input%ihwHO
           write(Input%chwHO,'(a4,i2)') 'hwHO',Input%ihwHO
           read(lin,'(A10,i2)')  ctemp,Input%iCOM
           if(Input%iCOM.eq.1) then
             if(.NOT. ALLOCATED(Input%ctpp)) ALLOCATE(character(len=7) :: Input%ctpp)
             Input%ctpp='_notpp_'
             if(lpr) write(*,'(a)') '     Center-Of-Mass Correction: 1B'
           elseif(Input%iCOM.eq.2) then
             if(.NOT. ALLOCATED(Input%ctpp)) ALLOCATE(character(len=5) :: Input%ctpp)
             Input%ctpp='_tpp_'
             if(lpr) write(*,'(a)') '     Center-Of-Mass Correction: 1B+2B'
           else
             if(lpr) write(*,*) '  Error: Center-of-Mass Correction is not set properly'
           endif
         endif

        if(Input%IntType.eq. 1 .and. Input%cIntID == 'GCN') &
     &  stop ' Interation Type is Wrong !' 
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'KB3') &
     &  stop ' Interation Type is Wrong !' 
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'USD') & 
     &  stop ' Interation Type is Wrong !' 

        read(lin,'(A10,I5)')   b(4),Input%IntIMSRG
        if(.NOT. ALLOCATED(Input%cFlow)) ALLOCATE(character(len=5) :: Input%cFlow)
        read(lin,'(A10,a5)')   b(5),Input%cFlow
        read(lin,'(A10,I5)')   b(5),Input%IntJT
        read(lin,'(A10,I5)')   b(6),Input%IsoSpin
   
        if(lpr) then
        if(Input%IntJT.eq.2)   print *,' -> Transform V2B from J-scheme to JTMT scheme '
        if(Input%IntJT.eq.1)   print *,' -> Decouple V2B from JTMT to m-scheme '
        if(Input%IntJT.eq.1)   print *,' -> Read V2B from  m-scheme '
        if(Input%IsoSpin.eq.0) print *,'    Allow the V2B breaking isospin symmetry '
        if(Input%IsoSpin.eq.1) print *,'    The V2B has isospin symmetry '
        endif


        read(lin,'(A10,I6)')   b(7),Input%nprot
        read(lin,'(A10,I6)')   b(8),Input%nneut
        read(lin,'(A10,I6)')   b(9),PNP%NFOM
        PNP%MPhi=PNP%NFOM             ! store the number of meshpoints for gauge angle 
        read(lin,'(A10,I6)')   b(10),AMP%NLEG_ALP
        read(lin,'(A10,I6)')   b(11),AMP%NLEG_BET
        read(lin,'(A10,I6)')   b(12),AMP%NLEG_GAM
! .............................

        if(PNP%NFOM .le. 1)   PNP%NFOM = 1
        if(AMP%NLEG_ALP*AMP%NLEG_BET*AMP%NLEG_GAM .ge. 1)   AMP%IsAMP = 1
        if(AMP%NLEG_ALP*AMP%NLEG_GAM.gt.1) then
           AMP%i3DAMP = 1
        else
           AMP%i3DAMP = 0
        endif

        read(lin,'(A10,i2)')    b(13),Input%iscale
        read(lin,'(A10,I20)')   ctemp,Input%ICons
        read(lin,'(A10,I20)')   b(14),Input%nq0i !Nuphi_ini
        read(lin,'(A10,I20)')   b(15),Input%nq0f !Nuphi_fin
        read(lin,'(A10,I20)')   b(16),Input%nq1i !nnnini
        read(lin,'(A10,I20)')   b(17),Input%nq1f !nnfin
        read(lin,'(A10,I6)')    b(18),Input%idens
        read(lin,'(A10,I6)')    b(19),Input%iRho3B
        read(lin,'(A10,a)')     ctemp,Input%vs4me3b
        if(lpr .and. Input%iRho3B .eq.1) &
       & print *, '    Computing Rho3B in Model Space:',Input%vs4me3b
        read(lin,'(A10,I6)')   b(21),Input%iE2
        read(lin,'(A10,I6)')   b(22),Input%icr

!  ............................
!        PNP%eps = 10.d0**(neps)
!       write(*,*) '   A truncation factor is adopted :',PNP%eps

        Input%Cscale = 10.d0**(Input%iscale)
        if(lpr .and. Input%iscale.ne.0)&
       &write(*,*) ' A scale factor is adopted for computing pfaffian:',Input%Cscale
       if(Input%nq0f.gt.nosmax .or. Input%nq1f.gt.nosmax) then
           if(lpr) print*,'increase nosmax in main file'
           stop
        end if

! ..............................

        ACONS_MV(1)= Input%nprot
        ACONS_MV(2)= Input%nneut
!   ................
        return
       end
