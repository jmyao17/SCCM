!______________________________________________________________________________
      subroutine readkernel_general(iq1,iq2,jprmi,jprma,jdf)
!..............................................................................
!     read matrix elements for Hill-Wheeler-Griffin solver to file lou      .
!..............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-z)

  401 format (10i5)
  402 format (8e15.8)

      print *, ' --> readkernel_general:'
!    ......................................................... 
      icheck = 0
      do jpr = jprmi,jprma,jdf
        if(jpr.eq.1) cycle
!       .................................... write matrix elements to tape 
        read (lou,401,end=403) jproj
        if(jproj.ne.jpr) stop 'Error: in reading kernels !'
        read (lou,402) bet1,gam1,bet2,gam2
!    .........................................................
      !  phi1 = Const%P00_mesh(iq1)
      !  phi2 = Const%P00_mesh(iq2)
!    .........................................................
        if(iv(jproj).gt.zero) then         ! for even spin
           k2max  =  jproj
        endif
        if(iv(jproj).lt.zero) then       ! for odd spin
           k2max  =  abs(jproj-1)
        endif
!    .........................................................
         do k1  =  -k2max, k2max, 1
         do k2  =  -k2max, k2max, 1
               read (lou,401) k1r,k2r
               !print *, jproj, k1r,k2r
               read (lou,402) Kernel%njkk(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%hjkk(iq1,iq2,jproj,k1,k2)

               read (lou,402) Kernel%nn(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%zz(iq1,iq2,jproj,k1,k2)

               read (lou,402) Kernel%J2(iq1,iq2,jproj,k1,k2)

               ! <J1=J2,K1 ||Q2|| J2,K2>
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,0),Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)

               read (lou,402) Kernel%q0p(iq1,iq2,jproj,k1,k2),Kernel%q0pc(iq1,iq2,jproj,k1,k2)
         enddo ! k2
         enddo ! k1                           

         do k2=  -k2max, k2max, 1
              do k1=  -k2max-1, k2max+1, 1
               ! <J1=J2+1,K1 ||Q2|| J2,K2>
                 read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,1),Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
              enddo ! k1
              do k1 =  -k2max-2, k2max+2, 1
               ! <J1=J2+2,K1 ||Q2|| J2,K2>
                 read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
              enddo ! k1
          enddo ! k2


         if(iq1.ne.iq2) then
         Kernel%q0p(iq2,iq1,jproj,k1,k2)     = Kernel%q0pc(iq1,iq2,jproj,k1,k2)
         Kernel%qpred(iq2,iq1,jproj,k1,k2,0) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
         do k2=  -k2max, k2max, 1
              do k1=  -k2max-1, k2max+1, 1
                 Kernel%qpred(iq2,iq1,jproj,k1,k2,1) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
              enddo ! k1
              do k1 =  -k2max-2, k2max+2, 1
                 Kernel%qpred(iq2,iq1,jproj,k1,k2,2) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
              enddo ! k1
          enddo ! k2
         endif

      icheck = icheck + 1
  99  enddo ! jproj
  403 continue
      return
      end 

!    ............................................
      subroutine table_general(jproj,iexst)
!    ............................................
!    iqk(iq,k)=1,2,.., NOQ*(2Jmax+1)
!    ............................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
      real*8 qtot(GCM%NOQ)

  201 format ('  J  K1 K2 iq     beta   gamma     qtot       P00       E  &
     &      ',           '   n^J(q,q)    <N>     <Z>     <J^2>')
  202 format (4i3,4f10.3,2f12.4,3f9.4)
      !write(*,201)
!     ................................. initialization
      maxmp = GCM%NOQ*(jproj*2+1)
      amas  = Nucl%nucleon(2) 
      r00   = r0*amas**third
      fac   = dsqrt(16*pi/5)
!     ................................ print diagonal (in q space)
!     matrix elements
      do iqk1 = 1,maxmp
           k1  = GCM%ik(iqk1,jproj)  
           iq1 = GCM%iq(iqk1,jproj)  
      do iqk2 = 1,maxmp
           k2  = GCM%ik(iqk2,jproj)   
           iq2 = GCM%iq(iqk2,jproj)
           ee  = GCM%hjkkqq(jproj,iqk1,iqk2)
           gg  = GCM%njkkqq(jproj,iqk1,iqk2)
           xnn = dreal(Kernel%nn(iq1,iq1,jproj,k1,k2))
           xj2 = dreal(Kernel%J2(iq1,iq1,jproj,k1,k2))
           xpp = dreal(Kernel%zz(iq1,iq1,jproj,k1,k2))
           qtot(iq1) = (3*amas*r00**2)/(4*pi)*Const%beta2t_mesh(iq1)*fac
           if(abs(gg).gt.CHOP) then
              GCM%Iexst(iqk1,iqk2) = log(abs(gg))
              if(k1.ne.k2 .or. iq1.ne.iq2) cycle
              write(*,202),jproj,k1,k2,iq1,Const%beta2t_mesh(iq1),  &
              Const%gamma2t_mesh(iq1),qtot(iq1),Const%P00_mesh(iq1),&
              ee/gg,gg,xnn,xpp,xj2
           endif
      enddo ! iq1
      enddo ! iq1

      return
      end

!    ............................................
      subroutine table(k1,k2,k1m,k2m,jproj,iexst)
!    ............................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
      real*8 qtot(GCM%NOQ)

  202 format (4i3,4f10.3,2f12.4,3f9.4)
!     ................................. initialization
      maxmp = GCM%NOQ
      amas  = Nucl%nucleon(2) 
      r00 = r0*amas**third
      fac = dsqrt(16*pi/5)
!     ................................ print diagonal (in q space) matrix elements
      do iq1 = 1,maxmp
           ee  = GCM%hjkkqq(jproj,iq1+k1m/2*maxmp,iq1+k2m/2*maxmp)
           gg  = GCM%njkkqq(jproj,iq1+k1m/2*maxmp,iq1+k2m/2*maxmp)
           xnn = dreal(Kernel%nn(iq1,iq1,jproj,k1,k2))
           xj2 = dreal(Kernel%J2(iq1,iq1,jproj,k1,k2))
           xpp = dreal(Kernel%zz(iq1,iq1,jproj,k1,k2))
           qtot(iq1) = (3*amas*r00**2)/(4*pi)*Const%beta2t_mesh(iq1)*fac
           if(abs(gg).gt.CHOP) &
           write(*,202),jproj,k1,k2,iq1,Const%beta2t_mesh(iq1),  &
                     Const%gamma2t_mesh(iq1),qtot(iq1),Const%P00_mesh(iq1),ee/gg,gg,xnn,xpp,xj2
      enddo ! iq1

!     ................................ print matrix of existing matrix elements
      do iq1 = 1,maxmp
      do iq2 = 1,maxmp
           iqk1=iq1+k1m/2*maxmp
           iqk2=iq2+k2m/2*maxmp
           gg  = GCM%njkkqq(jproj,iqk1,iqk2)
          if (abs(gg).gt.0.0d0) then
            GCM%Iexst(iqk1,iqk2) = log(abs(gg))
          endif
      enddo ! iq1
      enddo ! iq2

      return
      end
!______________________________________________________________________________
      subroutine readkernel_axial(iq1,iq2,jprmi,jprma,jdf)
!..............................................................................
!     read matrix elements for Hill-Wheeler-Griffin solver to file lou      .
!..............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-z)

  401 format (10i5)
  402 format (8e15.8)

!    ......................................................... 
      k1 = 0
      k2 = 0
      icheck = 0
      do jpr = jprmi,jprma,jdf
        if(jpr.eq.1) cycle
!       .................................... write matrix elements to tape 
        read (lou,401,end=403) jproj
        if(jproj.ne.jpr) stop 'Error: in reading kernels !'
        read (lou,402) bet1,gam1,bet2,gam2

        if(abs(Const%beta2t_mesh(iq1)-bet1).gt.0.01) stop ' beta1 is not correct !'
        if(abs(Const%gamma2t_mesh(iq1)-gam1).gt.1.0) stop ' gam1  is not correct !'
        if(abs(Const%beta2t_mesh(iq2)-bet2).gt.0.01) stop ' beta2 is not correct !'
        if(abs(Const%gamma2t_mesh(iq2)-gam2).gt.1.0) stop ' gam2  is not correct !'
!    .........................................................
        phi1 = Const%P00_mesh(iq1)
        phi2 = Const%P00_mesh(iq2)
!    .........................................................
!         if(abs(Kernel%njkk(jproj,k1,k2)).lt.1.d-8) cycle
               read (lou,401) k1r,k2r
               if(k1r.ne.0 .or. k2r.ne.0) &
               stop ' For axial states, only K=0 is considered'
               read (lou,402) Kernel%njkk(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%hjkk(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%nn(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%zz(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%J2(iq1,iq2,jproj,k1,k2)

!               if(abs(Input%nprot-Kernel%zz(iq1,iq2,jproj,k1,k2)).gt.0.1) &
!     &         write(*,'(4f6.2,a,2f6.2)') bet1,phi1,bet2,phi2, &
!     &         ' WARNING: <Z>=',Kernel%zz(iq1,iq2,jproj,k1,k2) 
!               if(abs(Input%nneut-Kernel%nn(iq1,iq2,jproj,k1,k2)).gt.0.1) &
!     &         write(*,'(4f6.2,a,2f6.2)') bet1,phi1,bet2,phi2, &
!     &         ' WARNING: <N>=',Kernel%nn(iq1,iq2,jproj,k1,k2) 

!               if((abs(bet1).gt.0.01 .and. abs(bet2).gt.0.01 ) .and. &
!     &             abs(jproj*(jproj+1)-Kernel%J2(iq1,iq2,jproj,k1,k2)).gt.0.1) &
!     &         write(*,'(4f6.2,a,i4,a,2f6.2)') bet1,phi1,bet2,phi2, &
!     &         ' WARNING: should be',jproj*(jproj+1),&
!     &                  ' but as ',Kernel%J2(iq1,iq2,jproj,k1,k2) 

               read (lou,402) Kernel%q0p(iq1,iq2,jproj,k1,k2),Kernel%q0pc(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,0),Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,1),Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
               !print *, jproj,k1,k2,Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
               !print *, jproj,k1,k2,Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
!     ....................................... qcp, qcpred were set to zero for iq1=iq2 cases in PKC code
!          q0p    = <J(f),K1,q1 ||e*r^2|| J(i),K2,q2>
!          qpred  = <J,J+1,J+2,K1,q1 ||e*Q_2|| J,K2,q2>
!          qcpred = <J,J+1,J+2,K1,q2 ||e*Q_2|| J,K2,q1>
!     .... exchange q1 and q2
          if(iq1.ne.iq2) then
            Kernel%qpred(iq2,iq1,jproj,k1,k2,0) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
            Kernel%qpred(iq2,iq1,jproj,k1,k2,1) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
            Kernel%qpred(iq2,iq1,jproj,k1,k2,2) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
            Kernel%q0p(iq2,iq1,jproj,k1,k2)   = Kernel%q0pc(iq1,iq2,jproj,k1,k2)
            Kernel%nn(iq2,iq1,jproj,k1,k2)    = Kernel%nn(iq1,iq2,jproj,k1,k2)
            Kernel%zz(iq2,iq1,jproj,k1,k2)    = Kernel%zz(iq1,iq2,jproj,k1,k2)
          endif
      icheck = icheck + 1
!       write(*,*) 'icheck=',icheck 
  99  enddo ! jproj
  403 continue

      if(icheck.eq.0) then
      write(12,*) elem, 'does not exist ...'
      do jproj = jprmi,jprma,jdf
         do k1  =  0, jproj, 2
         do k2  =  0, jproj, 2
            Kernel%njkk(iq1,iq2,proj,k1,k2)   = 0.d0
            Kernel%hjkk(iq1,iq2,jproj,k1,k2)  = 0.d0
            Kernel%qpred(iq2,iq1,jproj,k1,k2,0:2) = 0.d0
            Kernel%qpred(iq1,iq2,jproj,k1,k2,0:2) = 0.d0
            Kernel%q0p(iq2,iq1,jproj,k1,k2)   = 0.d0
            Kernel%q0p(iq1,iq2,jproj,k1,k2)   = 0.d0
            Kernel%nn(iq1,iq2,jproj,k1,k2)= 0.d0
            Kernel%nn(iq2,iq1,jproj,k1,k2)= 0.d0
            Kernel%zz(iq2,iq1,jproj,k1,k2)= 0.d0
            Kernel%zz(iq1,iq2,jproj,k1,k2)= 0.d0
          enddo ! k2 
          enddo ! k1 
      enddo ! jproj
      endif
      return
      end 



!.......................................................................
      subroutine ReadCoordinators()
      use vaphfb_par
      implicit none
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
         if(.NOT. ALLOCATED(Const%hw_mesh))       ALLOCATE(Const%hw_mesh(1:NQ))
         do iq=1,NQ
            read(15,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),&
            Const%P00_mesh(iq),Const%hw_mesh(iq)
            write(*,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),Const%P00_mesh(iq)
            call Filename4wfs(Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),Const%P00_mesh(iq),Const%hw_mesh(iq))
            File%cwf(iq) = File%wf    ! name of files for wfs
         enddo !
          Input%NGCM=NQ
!         print *,'Input%NGCM=',Input%NGCM
         if(Input%nq0f.gt.NQ .or. Input%nq1f.gt.NQ) &
     &   stop ' Error: Input%nq0f or Input%nq1f are larger than the total &
         &number of configurations!'
      end

   
!
