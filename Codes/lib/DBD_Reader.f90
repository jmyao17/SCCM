      subroutine Reader(lpr)
!    ...................
!    eta1     -> eta
!    Nitermax - > itermax
!    ...................

       USE VAPHFB_Par
       implicit none 
       logical lpr
       character*10 ctemp,a(50),b(50)
       integer ii,isum,IsGCM

        open(5,file='input.dat',status='old')
!       ..................READ THE INPUT PARAMETERS
        read(5,'(A10,i1)')    b(1),Input%IntType
        write(6,'(A10,i1)')    b(1),Input%IntType
!       ................. shell model int
        if(Input%IntType.eq.0) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20):: Input%cIntID)
           read(5,'(A10,a18)')  b(2),Input%cIntID
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6):: Input%cValID)
           read(5,'(A10,a5)')   b(3),Input%cValID

            write(*,'(a,a18)') 'Interaction Type: ',Input%cIntID
            write(*,'(a,a5)') 'Model Space: ',trim(Input%cValID)
         endif

!       ................. wfs by chiral NN int
         if(Input%IntType.eq.1) then
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6)  :: Input%cValID)
           read(5,'(A10,a20)')  b(2),Input%cIntID
           read(5,'(A10,a6)')   b(3),Input%cValID
         endif
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'GCN') &
     &  stop ' Interation Type is Wrong !' 
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'KB3') &
     &  stop ' Interation Type is Wrong !' 
        if(Input%IntType.eq. 1 .and. Input%cIntID == 'USD') & 
     &  stop ' Interation Type is Wrong !' 

! .............................

        read(5,'(A10,I6)')   b(4),Input%NME_Type
        read(5,'(A10,i2)')  ctemp,Input%ihwHO
        print *, 'DBD_TBME: ', Input%NME_Type
        print *, 'hw: ', Input%ihwHO

        if(.NOT. ALLOCATED(Input%cFlow)) ALLOCATE(character(len=6) :: Input%cFlow)
        read(5,'(A10,a6)')   ctemp,Input%cFlow
        read(5,'(A10,i1)')   ctemp,Input%IsHFB

        read(5,'(A10,I6)')   b(5),Input%nprot
        read(5,'(A10,I6)')   b(6),Input%nneut

        ACONS_MV(1)= Input%nprot
        ACONS_MV(2)= Input%nneut
        print *, 'Proton: ',  Input%nprot
        print *, 'Neutron: ', Input%nneut

        read(5,'(A10,I6)')    b(7),Input%idens
        read(5,'(A10,I6)')    b(8),PNP%NFOM
        PNP%MPhi = PNP%NFOM  
        read(5,'(A10,I6)')    b(9),AMP%NLEG_ALP
        read(5,'(A10,I6)')    ctemp,AMP%NLEG_BET
        read(5,'(A10,I6)')    ctemp,AMP%NLEG_GAM
        if(PNP%NFOM .le. 1)   PNP%NFOM = 1
        if(AMP%NLEG_ALP*AMP%NLEG_BET*AMP%NLEG_GAM .ge. 1)   AMP%IsAMP = 1

        read(5,'(A10,I6)')    a(15),Input%iTD_Calc
        read(5,'(A10,I6)')    ctemp,IsGCM

        if (IsGCM.eq.1) then
          print *, ' calculate transition density with config. mixing'
          return
        endif
        print *, ' Computing two-body transition density'

        read(5,'(A10,i6)')    a(10),Input%iscale
        read(5,'(A10,I20)')   b(14),Input%nqas !Nuphi_ini
        read(5,'(A10,I20)')   b(15),Input%nqae !Nuphi_fin
        read(5,'(A10,I20)')   b(16),Input%nqbs !nnnini
        read(5,'(A10,I20)')   b(17),Input%nqbe !nnfin


        if(lpr) write(*,*) Input%iscale,Input%iTD_Calc
        Input%Cscale = 10.d0**(Input%iscale)
        if(lpr) &
        write(*,*) ' A scale factor is adopted for computing pfaffian:',&
                   Input%Cscale

        if(Input%Nuphi_ini.gt.Nuphimax .or. Input%Nuphi_fin.gt.Nuphimax) then
          print*,'increase Nuphimax in main file'
          stop
        end if

!......... printout

       if(lpr) then
          write(6,'(A)')   '....................' 
          write(6,'(A)')   ' Input Parameters   ' 
          write(6,'(A)')   '....................' 
       endif
! ..............................
        return
       end
