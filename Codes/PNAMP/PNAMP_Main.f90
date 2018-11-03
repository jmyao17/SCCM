
        program PNAMP_MPI
 
        ! include mpi module
        use mpi
        use vaphfb_par
        implicit none
        logical lpr
        integer Nphi1,Nphi2,Nphi2s,ic,indice
        ! Data declaration for mpi
        integer :: ierr
        integer :: rank 
        integer :: nprocs 

       ! initialize mpi
       call MPI_INIT(ierr)
       ! setup communicator size
       call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
       ! setup ranks/IDs for each process
       call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

       !print *, "hello from",rank
       if(rank.eq.0) lpr =.true. 
       if(rank.ne.0) lpr =.false. 


       call PNAMP_Default() 
       call Reader(lpr)

       call Prep(lpr)
       call ReadCoordinators(lpr)
       print *, 'model space'
       call Model_space()
       print *, 'const'
       call ZCONST
       print *, 'gfv'
       call gfv
       print *, 'ME2B_Base'
       call ME2B_Base_Full()

       if(Input%iRho3B.eq.1) then
        call Rho3B_1B_Model_Space(lpr)
        call Rho3B_2B_Model_Space(lpr)
        call ME3B_Base_Full(lpr)
       endif
       call Hamiltonian(rank,lpr)
       call Isospin_t3()
       call SingleParticleME()
       !call Reduced_ME1B()
       call Gaussmesh()
!     --------------------------------------
        print *, ' -> Information about the output files: '
        if(lpr) then
         write(311,*)Input%nq0i,Input%nq0f
         write(311,*)Input%nq1i,Input%nq1f
        endif


!        determine the configurations to be calculated
         ic=0
         DO Nphi1 =Input%nq0i,Input%nq0f       ! <phi1| ... |phi2>
             if(Input%ICons.eq.1) then
                Nphi2s=Input%nq1i
             else                 ! using symmetry, only half is calculated 
                Nphi2s=Nphi1
             endif
         DO Nphi2 =Nphi2s,Input%nq1f       ! <phi1| ... |phi2>
            ic=ic+1
            Kernel%Nphi1(ic)=Nphi1
            Kernel%Nphi2(ic)=Nphi2
         enddo
         enddo
         Kernel%Nelem    =ic
         if(rank.eq.0) then 
            print *, 'Total No. of Kernels to be calculated:',Kernel%Nelem
            print *, 'Total No. of threads/ranks:',nprocs
         endif
!            Num_Loops = Kernel%Nelem/nprocs+1  
         do ic=1,Kernel%Nelem 
           ! the process with the rank runs the calculation for the
           ! corresponding
           ! determine loops

           if( rank .eq. mod(ic,nprocs)) then

               Nphi1    = Kernel%Nphi1(ic)
               Nphi2    = Kernel%Nphi2(ic)
               PNP%NFOM = PNP%MPhi 
 
               
               print *, 'Rank:',rank,' is running ...'
               !print *, 'Nphi1=',Nphi1,' Nphi2=',Nphi2
               !print *, ' The ', ic,'/',Kernel%Nelem, &
               !' elements is in commputation'

!  ..............................................
!    kernel
!  ..............................................

        call Filename4Kernels(AMP%i3DAMP,Nphi1,Nphi2)
        if(AMP%i3DAMP.eq.0)  then   ! 1DAMP for axial case
            write(*,10) rank,Nphi1,Nphi2,File%elem
            !write(*,*) ' Density matrix elements are written to '
            !write(*,*) File%Rho1B,File%Rho2B,File%Rho3B
            call Kernel4Axial(lpr,Nphi1,Nphi2,HO%NLEV)
            !call Kernel4general(.true.,Nphi1,Nphi2,HO%NLEV)
            open(lou,file=File%elem,status='unknown')
            call post_axial(0,Nphi1,Nphi2,0,JJmax,2)
        else                         ! 3DAMP for non-axial cases
            print *, ' ---> 3DAMP ...'
            write(*,*) ' Kernels are written to '
            write(*,*) File%elem
            if(Input%icr.eq.0) call Kernel4Triaxial(.true.,Nphi1,Nphi2,HO%NLEV)
            if(Input%icr.eq.1) call Kernel4Crankingx(.true.,Nphi1,Nphi2,HO%NLEV)
            if(Input%icr.eq.2) call Kernel4general(.true.,Nphi1,Nphi2,HO%NLEV)
            open(lou,file=File%elem,status='unknown')
            call post_general(1,Nphi1,Nphi2,0,JJmax,1)
        endif     

        endif  ! end rank
9       END DO !end do loop ic 
10    format('rank:',i4,'<q',i3,'|O|',i3,'> to be saved into ',a) 
!   ............................. 
        ! finalize MPI
        print *, 'Rank:',rank,' is finished ...'
        call MPI_FINALIZE(ierr)

99     stop  
       end
