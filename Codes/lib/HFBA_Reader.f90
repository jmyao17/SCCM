     subroutine Reader(lpr)
!    ...................
!    READ THE INPUT PARAMETERS
!    eta1     -> eta
!    Nitermax - > itermax
!    ...................

       USE VAPHFB_Par
       implicit none 
       logical lpr
       character*10 ctemp,b(0:50)
       character*30 ctemp30
       integer itemp,ii,it,isum

       if(lpr) then
       print *, ' ======================== '
       print *, ' Features of Calculations '
       print *, ' ======================== '
       endif

        open(lin,file='hfb.dat',status='old')
        read(lin,'(A10,i1)')    b(0),Input%NPMix
        read(lin,'(A10,i1)')    ctemp,Input%IPMix
        read(lin,'(A10,i1)')    ctemp,Input%KMix
        read(lin,'(A10,i1)')    ctemp,Input%IsHFB
        read(lin,'(A10,i1)')    b(1),Input%IntType
        do it=0,1
           if(Input%IsHFB.eq.0) HFB%ide(it) = 0
           if(Input%IsHFB.eq.1) HFB%ide(it) = 1
        enddo !it

        if(lpr) then
        if(Input%IsHFB.eq.0) then
           print *, ' -> Hartree-Fock Calculation'
        else
           print *, ' -> Hartree-Fock-Bogoliubov Calculation '
        endif

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
!       ................. shell model int
        if(Input%IntType.eq.0) then
           read(lin,'(A10,i1)')    ctemp,Input%InME3B
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6) :: Input%cValID)
           read(lin,'(A10,a18)')  b(2),Input%cIntID
           read(lin,'(A10,a5)')   b(3),Input%cValID
           read(lin,'(A10,i2)')  ctemp,Input%ihwHO
           read(lin,'(A10,i2)')  ctemp,Input%iCOM
           write(Input%chwHO,'(a4,i2)') 'hwHO',Input%ihwHO 
!        ..........................
           if(lpr) write(*,'(a,a5)') '     Model Space: ',trim(Input%cValID)
           if(lpr) write(*,'(a,a18)') '     Interaction Type: ',Input%cIntID
         endif
!       ................. chiral interation
         if(Input%IntType.eq.1) then
           read(lin,'(A10,i1)')    ctemp,Input%InME3B
           if(.NOT. ALLOCATED(Input%cIntID)) ALLOCATE(character(len=20) :: Input%cIntID)
           if(.NOT. ALLOCATED(Input%cValID)) ALLOCATE(character(len=6)  :: Input%cValID)
           read(lin,'(A10,a20)')  b(2),Input%cIntID
           read(lin,'(A10,a6)')   b(3),Input%cValID
           read(lin,'(A10,i2)')  ctemp,Input%ihwHO
           read(lin,'(A10,i2)')  ctemp,Input%iCOM
           write(Input%chwHO,'(a4,i2)') 'hwHO',Input%ihwHO 
           if(lpr) then
           write(*,'(a,a20)') '     Interaction Type: ',Input%cIntID
           write(*,'(a,a6)') '     Model Space: ',Input%cValID
           write(*,'(a,a6)') '     Frequency of HO: ',Input%chwHO
           endif
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

         if(Input%IntType.eq. 1 .and. Input%cIntID == 'GCN') stop ' Interation Type is Wrong !' 
         if(Input%IntType.eq. 1 .and. Input%cIntID == 'KB3') stop ' Interation Type is Wrong !' 
         if(Input%IntType.eq. 1 .and. Input%cIntID == 'USD') stop ' Interation Type is Wrong !' 
        

         if(.NOT. ALLOCATED(Input%cFlow)) ALLOCATE(character(len=5) :: Input%cFlow)
        read(lin,'(A10,a5)')   b(5),Input%cFlow
        read(lin,'(A10,I5)')    b(5),Input%IntJT
        read(lin,'(A10,I5)')    b(6),Input%IntIMSRG
        read(lin,'(A10,I5)')    b(7),Input%IsoSpin
        read(lin,'(A10,I6)')    b(8),Input%nprot
        read(lin,'(A10,I6)')    b(9),Input%nneut
        print *, Input%nprot,Input%nneut
        read(lin,'(A10,I5)')    b(10),Input%itermax
        read(lin,'(A10,a5)')    ctemp,Opt%method   ! optimize method 
        read(lin,'(A10,F20.6)') b(11),Opt%lr       ! learning rate or step size
        print *, ' -> Optimizer adopted  '
        write(*,*) '     ',Opt%method
        write(*,'(a,f10.5)') '     Learning rate:',Opt%lr
        if(Opt%method =='SDM') &
        write(*,'(a,f10.5)') '     Momentum Para:',Opt%mp
        read(lin,'(A10,F20.6)') b(12),Input%tolcons
        read(lin,'(A10,F20.6)') b(13),Input%tolgrad
        read(lin,'(A10,I5)')    b(14),Input%inwf
        read(lin,'(A10,I5)')    b(15),PNP%NFOM

        if(PNP%NFOM .le. 1)   PNP%NFOM = 1

!    ............ printout
      if(lpr) then
        if(Input%IntJT.eq.2)   print *,' -> Transform V2B from J-scheme to JTMT scheme '
        if(Input%IntJT.eq.1)   print *,' -> Decouple V2B from JTMT to m-scheme '
        if(Input%IntJT.eq.1)   print *,' -> Read V2B from  m-scheme '
        if(Input%IsoSpin.eq.0) print *,'    Allow the V2B breaking isospin symmetry '
        if(Input%IsoSpin.eq.1) print *,'    The V2B has isospin symmetry '

        print *,' -> Information about the Nucleus: '
        write(*,'(a,i3)') '     Proton: ',Input%nprot
        write(*,'(a,i3)') '     Neutron:',Input%nneut

        print *,' -> Information about the Iteration: '
        write(*,'(a,i3)')   '     Max. Iteration:',Input%itermax
        write(*,'(a,f8.5)') '     Step Size of Gradient Method:',Opt%lr
        write(*,'(a,f8.5)') '     Required precision:',Input%tolgrad
        if(PNP%NFOM.gt.1) print *,' -> PNP-VAPHFB method '
        if(PNP%NFOM.le.1) print *,' -> HFB method '
      endif
!   ............. for constraint calculation
        read(lin,*)
        read(lin,'(A10,f8.5)') ctemp,Input%bkick
        read(lin,'(A10,I5)')   b(16),Const%iQB
        read(lin,'(A10,I5)')   b(17),Const%NCONS


        isum=0
        if(lpr) write(6,*) '.......................'
        if(lpr) print *,' -> Constraint Parameters '
        do ii=3,NCONSMAX  ! 18 
            read(lin,'(A10,I1,F20.6)') &
     &           b(15+ii),Icons_y_n(ii),ACONS_MV_aux(ii)
            if(Icons_y_n(ii).eq.1) then
              if(lpr) write(*,'(a,a10,a,f10.5)') '     ', b(15+ii),' =>',ACONS_MV_aux(ii) 
            endif
             isum = isum + Icons_y_n(ii)
        end do

        if(lpr) write(6,*) '.......................'
        if(isum.ne.Const%NCONS) then
           write(6,*) 'The sum of constraints does not correspond &
           & with NCONS, fix it in the input file'
           stop
        end if

        ACONS_MV(1)= Input%nprot
        ACONS_MV(2)= Input%nneut
!       ........ plus contraints on N and Z
        Const%NCONS=COnst%NCONS+2 !adding proton/neutron constraint

!   ................
        return
       end
