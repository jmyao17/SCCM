!______________________________________________________________________________
      subroutine post_axial(icase,iq1,iq2,jprmi,jprma,jdf)

!..............................................................................
!     write matrix elements for Hill-Wheeler-Griffin solver to file lou
!..............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-z)

  401 format (10i5)
  402 format (8e15.8)

      do jproj = jprmi,jprma,jdf
         if(jproj.eq.1) cycle 
          write (lou,401) jproj
          write (lou,402) Const%beta2t_mesh(iq1),Const%gamma2t_mesh(iq1),&
     &                    Const%beta2t_mesh(iq2),Const%gamma2t_mesh(iq2)
!    .........................................................
          k1 = 0
          k2 = 0
          write (lou,401) k1,k2
          write (lou,402) Kernel%njkk(jproj,k1,k2), &
     &                    Kernel%hjkk(jproj,k1,k2)
          write (lou,402) Kernel%nn(jproj,k1,k2), &
     &                    Kernel%zz(jproj,k1,k2)
! .................... newly added for J^2
          write (lou,402) Kernel%J2(jproj,k1,k2)
          write (lou,402) Kernel%q0p(jproj,k1,k2),Kernel%q0pc(jproj,k1,k2)
          write (lou,402) Kernel%qpred(jproj,k1,k2,0),Kernel%qcpred(jproj,k1,k2,0)
          write (lou,402) Kernel%qpred(jproj,k1,k2,1),Kernel%qcpred(jproj,k1,k2,1)
          write (lou,402) Kernel%qpred(jproj,k1,k2,2),Kernel%qcpred(jproj,k1,k2,2)
 
!     ......................................
  99  enddo ! jproj
      return
      end 

!______________________________________________________________________________
      subroutine post_general(icase,iq1,iq2,jprmi,jprma,jdf)

!.............................................................................
!     write matrix elements for Hill-Wheeler-Griffin solver to file lou      .
!.............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-z)

  401 format (10i5)
  402 format (8e15.8)

      do jproj = jprmi,jprma,jdf
         if(jproj.eq.1) cycle 
!       .................................... write matrix elements to tape 
        write (lou,401) jproj
        write (lou,402) Const%beta2t_mesh(iq1),Const%gamma2t_mesh(iq1), &
     &                  Const%beta2t_mesh(iq2),Const%gamma2t_mesh(iq2)
!    .........................................................
        if(iv(jproj).gt.zero) then         ! for even spin
           k2max  =  jproj
        endif
        if(iv(jproj).lt.zero) then       ! for odd spin
           k2max  =  abs(jproj-1)
        endif
!    .........................................................
         do k2  =  -k2max, k2max, 1
         do k1  =  -k2max, k2max, 1
               write (lou,401) k1,k2
               write (lou,402) Kernel%njkk(jproj,k1,k2),                    &
     &                         Kernel%hjkk(jproj,k1,k2) !/Kernel%njkk(jproj,k1,k2)
               write (lou,402) Kernel%nn(jproj,k1,k2), &
     &                         Kernel%zz(jproj,k1,k2)
! .................... newly added for J^2
               write (lou,402) Kernel%J2(jproj,k1,k2)
               write (lou,402) Kernel%qpred(jproj,k1,k2,0),Kernel%qcpred(jproj,k1,k2,0)
               write (lou,402) Kernel%q0p(jproj,k1,k2),Kernel%q0pc(jproj,k1,k2)
              ! added for test 
             ! write (lou,402) Kernel%qpred(jproj,k1,k2,2),Kernel%qcpred(jproj,k1,k2,2)
              !write (*,'(i4,4f12.8)') jproj,Kernel%qpred(jproj,k1,k2,0),Kernel%qcpred(jproj,k1,k2,0)
          enddo ! k2
         enddo ! k1

         !  <J2,K2|| Q2 || J1,K1>, where J2=J1,J1+1,J1+2
         do k2  =  -k2max, k2max, 1
              do k1  =  -k2max-1, k2max+1, 1
                 write (lou,402) Kernel%qpred(jproj,k1,k2,1),Kernel%qcpred(jproj,k1,k2,1)
              enddo ! k1
              do k1  =  -k2max-2, k2max+2, 1
                 write (lou,402) Kernel%qpred(jproj,k1,k2,2),Kernel%qcpred(jproj,k1,k2,2)
              enddo ! k1
          enddo ! k2

!     ......................................
  99  enddo ! jproj
      return
      end 

!.......................................................................

      subroutine ReadCoordinators(lpr)
      use vaphfb_par
      implicit none
      logical lpr
      integer NQ,iq
!    .......... read coordinates
         open(15,file='betgam.dat',status='old')
         read(15,*)
         read(15,*)
         read(15,*)
         read(15,*)
!     ............
         read(15,*) NQ
         if(.NOT. ALLOCATED(Const%beta2t_mesh))   ALLOCATE(Const%beta2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%gamma2t_mesh))  ALLOCATE(Const%gamma2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%P00_mesh))      ALLOCATE(Const%P00_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%hw_mesh))      ALLOCATE(Const%hw_mesh(1:NQ))
         do iq=1,NQ
            read(15,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),&
                       Const%P00_mesh(iq), Const%hw_mesh(iq)
            call Filename4wfs(Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),&
                              Const%P00_mesh(iq), Const%hw_mesh(iq))
            File%cwf(iq) = File%wf    ! name of files for wfs
            if(lpr) print *,File%cwf(iq)
         enddo !
          Input%NGCM=NQ
            if(lpr) print *,'Input%NGCM=',Input%NGCM
         if(Input%nq0f.gt.NQ .or. Input%nq1f.gt.NQ) &
     &   stop ' Error: Input%nq0f or Input%nq1f are larger than the total &
         &number of configurations!'
      end

!________________________________________________________________________



