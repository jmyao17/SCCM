        subroutine ReadGeneCoord()

        use VAPHFB_PAR
        implicit none
        integer NQ,iq
!     ........................................
!     constraints on a set of states with diff. beta,gamma,p00
!     ........................................
      if(Const%iQB.eq.1)  then
         open(5,file='betgam.dat',status='old')
         read(5,*)
         read(5,*)
         read(5,*)
         read(5,*)
!     ............
         read(5,*) Const%NQ
         NQ=Const%NQ
         if(NQ.lt.1) stop 'Error: Number of Configs is set improperly'
         if(.NOT. ALLOCATED(Const%beta2t_mesh))   ALLOCATE(Const%beta2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%gamma2t_mesh))  ALLOCATE(Const%gamma2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%P00_mesh))      ALLOCATE(Const%P00_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Q20t_mesh))     ALLOCATE(Const%Q20t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Q22t_mesh))     ALLOCATE(Const%Q22t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Etot))          ALLOCATE(Const%Etot(1:NQ))
         if(.NOT. ALLOCATED(Const%hw))          ALLOCATE(Const%hw(1:NQ))
         do iq=1,NQ
            read(5,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq), &
     &                Const%P00_mesh(iq), &
     &                Const%hw(iq)
            Const%Q20t_mesh(iq) = Const%beta2t_mesh(iq) &
     &                           *dcos(Const%gamma2t_mesh(iq)*pi/180.0)/Const%Q2BA
            Const%Q22t_mesh(iq) = Const%beta2t_mesh(iq) &
     &            *dsin(Const%gamma2t_mesh(iq)*pi/180.0)/(Const%Q2BA*sqrt(2.d0))

!           print*,'Q20t_mesh',Const%Q20t_mesh(iq), 'Q22t_mesh=',Const%Q22t_mesh(iq)
!           print*, 'coeff=', Const%Q2BA
         enddo !
!     ......................................... read Q20t,Q22t
       else
         NQ = 1
         if(.NOT. ALLOCATED(Const%beta2t_mesh))   ALLOCATE(Const%beta2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%gamma2t_mesh))  ALLOCATE(Const%gamma2t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Q20t_mesh))     ALLOCATE(Const%Q20t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Q22t_mesh))     ALLOCATE(Const%Q22t_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%P00_mesh))   ALLOCATE(Const%P00_mesh(1:NQ))
         if(.NOT. ALLOCATED(Const%Etot))          ALLOCATE(Const%Etot(1:NQ))
         Const%Q20t_mesh(1)  = acons_mv_aux(3)
         Const%Q22t_mesh(1)  = acons_mv_aux(4)
         Const%P00_mesh(1)   = acons_mv_aux(12)

         Const%beta2t_mesh(1)  = (acons_mv_aux(3)**2 + ACONS_MV_aux(4)**2/2.d0)*Const%Q2BA
         if(abs(ACONS_MV_aux(3)).gt.1.d-8) then
           Const%gamma2t_mesh(1) = 180.d0*atan(sqrt(2.d0)*acons_mv_aux(4)/ACONS_MV_aux(3))/pi
         else
           Const%gamma2t_mesh(1) = 0.d0
         endif
      endif

         return
         end


       
       subroutine SetConstraintValue(iq,NLEV)
       USE VAPHFB_PAR
       implicit none
       integer iq4_const,ii,iq,NLEV


       ipair_const=0     ! counting number of # of pairing constraint terms
       iq2_const=0        ! counting number of # of quadrupole constraint terms
       iq4_const=0        ! counting number of # of hexapole constraint terms
       ind_con=2         ! counting number of total # of constraint terms
!      ................. constraints on Iso-scalar Q2
       if (Icons_y_n(3).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q20t(ii)
          end do
          acons_mv(ind_con+1)=Const%Q20t_mesh(iq)
          ind_con  = ind_con+1
          iq2_const= iq2_const + 1

!     .......................... stored the constraint on total mass Q20
        qtc20 = Const%Q20t_mesh(iq)
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(4).eq.1) then
           do ii=1,NLEV**2
              ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q22t(ii)
           end do
           acons_mv(ind_con+1)=Const%Q22t_mesh(iq)
           ind_con  = ind_con+1
           iq2_const= iq2_const + 1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(5).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q21t(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(5)
          ind_con=ind_con+1
          iq2_const= iq2_const + 1
       end if

!      ................. constraints on Iso-vector Q2
       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(6).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q20m(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(6)
          ind_con=ind_con+1
          iq2_const= iq2_const + 1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(7).eq.1) then
           do ii=1,NLEV**2
              ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q22m(ii)
           end do
           acons_mv(ind_con+1)=acons_mv_aux(7)
           ind_con=ind_con+1
           iq2_const= iq2_const + 1
       end if
       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(8).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q21m(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(8)
           ind_con=ind_con+1
           iq2_const= iq2_const + 1
      end if

!   ................ constraints on Jx,Jy,Jz

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(9).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%AJX_ME(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(9)
          print *, 'constraint on <Jx>:',acons_mv_aux(9)
          ind_con=ind_con+1
       end if
       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(10).eq.1) then
           do ii=1,NLEV**2
              ACONS_ME(ii+ind_con*NLEV**2)=cME1B%AJY_ME(ii)
           end do
           acons_mv(ind_con+1)=acons_mv_aux(10)
           ind_con=ind_con+1
       end if
       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(11).eq.1) then
           do ii=1,NLEV**2
              ACONS_ME(ii+ind_con*NLEV**2)=cME1B%AJZ_ME(ii)
           end do
           acons_mv(ind_con+1)=acons_mv_aux(11)
           ind_con=ind_con+1
       end if


!     .............. constraints on different types of pairing correlations

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(12).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_00_10_me(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(12)
         ind_con=ind_con+1
         ipair_const=ipair_const+1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(13).eq.1) then
           do ii=1,NLEV**2
              ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_1m1_00_me(ii)
           end do
           acons_mv(ind_con+1)=acons_mv_aux(13)
           ind_con=ind_con+1
           ipair_const=ipair_const+1
        end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(14).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_1p1_00_me(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(14)
          ind_con=ind_con+1
          ipair_const=ipair_const+1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(15).eq.1) then
         do ii=1,NLEV**2
            ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_10_00_me(ii)  ! P_TMT_JMJ_me
         end do
         acons_mv(ind_con+1)=acons_mv_aux(15)
         ind_con=ind_con+1
         ipair_const=ipair_const+1
       end if 

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(16).eq.1) then
         do ii=1,NLEV**2
            ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_00_1m1_me(ii)
         end do
         acons_mv(ind_con+1)=acons_mv_aux(16)
         ind_con=ind_con+1
         ipair_const=ipair_const+1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(17).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%P_00_1p1_me(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(17)
          ind_con=ind_con+1
          ipair_const=ipair_const+1
       end if

       if(ind_con.eq.NCONSMAX) goto 10
!    .... hexapole
       if (Icons_y_n(18).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q40t(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(18)
          ind_con=ind_con+1
          iq4_const=iq4_const+1
       end if
       if(ind_con.eq.NCONSMAX) goto 10
       if (Icons_y_n(19).eq.1) then
          do ii=1,NLEV**2
             ACONS_ME(ii+ind_con*NLEV**2)=cME1B%Q40m(ii)
          end do
          acons_mv(ind_con+1)=acons_mv_aux(19)
          ind_con=ind_con+1
          iq4_const=iq4_const+1
       end if
  10  continue
       iJ_const = ind_con-ipair_const-iq2_const-iq4_const-2 
       write(*,'(a44)')    '   Number of constraints on N&Z    :  2'
       write(*,'(a40,i3)') '   Number of constraints on Q2     :',iq2_const
       write(*,'(a40,i3)') '   Number of constraints on Q4     :',iq4_const
       write(*,'(a40,i3)') '   Number of constraints on Pairing:',ipair_const
       write(*,'(a40,i3)') '   Number of constraints on Ji     :',iJ_const
       write(*,'(a40,i3)') '   Total number of constraints     :',ind_con

       return
       end

!  ................................
        subroutine ME1B4Constraint()
!  ................................
        USE VAPHFB_PAR
        implicit none
        integer ii,i,j
!   ......................... allocate memory

      if(.NOT. ALLOCATED(cME1B%r2))  ALLOCATE(cME1B%r2(1:HO%NLEV*HO%NLEV))

      if(.NOT. ALLOCATED(cME1B%Q2_2t))  ALLOCATE(cME1B%Q2_2t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q2_1t))  ALLOCATE(cME1B%Q2_1t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q20t))  ALLOCATE(cME1B%Q20t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q21t))  ALLOCATE(cME1B%Q21t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q22t))  ALLOCATE(cME1B%Q22t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q2_2m))  ALLOCATE(cME1B%Q2_2m(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q2_1m))  ALLOCATE(cME1B%Q2_1m(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q20m))  ALLOCATE(cME1B%Q20m(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q21m))  ALLOCATE(cME1B%Q21m(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q22m))  ALLOCATE(cME1B%Q22m(1:HO%NLEV*HO%NLEV))

      if(.NOT. ALLOCATED(cME1B%Q30t))  ALLOCATE(cME1B%Q30t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q40t))  ALLOCATE(cME1B%Q40t(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(cME1B%Q40m))  ALLOCATE(cME1B%Q40m(1:HO%NLEV*HO%NLEV))

!  ........... radii
       call Rsq_ME(cME1B%r2)

       if(Const%iQBType.eq.0) then
!  ................. Kumar-Bagnar definition
         call q2mume_KB(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)
         call q2mumeIso_KB(cME1B%Q2_2m,cME1B%Q2_1m,cME1B%Q20m,cME1B%Q21m,cME1B%Q22m)
         call q3mume_KB(cME1B%Q30t)
!  ................. standard defintion
       else
!          call q2mume(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)
          call q2mume_ME1B(cME1B%Q2_2t,cME1B%Q2_1t,cME1B%Q20t,cME1B%Q21t,cME1B%Q22t)
          call q2mume_IsoV(cME1B%Q2_2m,cME1B%Q2_1m,cME1B%Q20m,cME1B%Q21m,cME1B%Q22m)
          call q3mume(cME1B%Q30t)
          call q40me_ME1B(cME1B%Q40t,0)
          call q40me_ME1B(cME1B%Q40m,1)
       endif

!       ..........angular momentum operators
!        write(*,*) 'Main: JZ_JZ2'

        if(.NOT. ALLOCATED(cME1B%AJX_ME))  ALLOCATE(cME1B%AJX_ME(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%AJY_ME))  ALLOCATE(cME1B%AJY_ME(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%AJZ_ME))  ALLOCATE(cME1B%AJZ_ME(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%AJX2_ME))  ALLOCATE(cME1B%AJX2_ME(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%AJY2_ME))  ALLOCATE(cME1B%AJY2_ME(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%AJZ2_ME))  ALLOCATE(cME1B%AJZ2_ME(1:HO%NLEV*HO%NLEV))

        call JZ_JZ2_ME(cME1B%AJX_ME,cME1B%AJY_ME,cME1B%AJZ_ME,cME1B%AJX2_ME,cME1B%AJY2_ME,cME1B%AJZ2_ME,HO%NLEV)

!        do i=1,HO%NLEV
!        do j=1,HO%NLEV
!           if(abs(cME1B%AJX2_ME(i+j*(HO%NLEV-1))).gt.1.d-4) &
!      &     write(700,'(2i3,f10.5)') i,j,cME1B%AJX2_ME(i+j*(HO%NLEV-1))
!        enddo
!        enddo

        !pairs
!       ..... iso-scalar pairing: np (T=0; J=1)
        if(.NOT. ALLOCATED(cME1B%P_00_1m1_me))  ALLOCATE(cME1B%P_00_1m1_me(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%P_00_10_me))  ALLOCATE(cME1B%P_00_10_me(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%P_00_1p1_me))  ALLOCATE(cME1B%P_00_1p1_me(1:HO%NLEV*HO%NLEV))

        call Pairs_TMT_JMJ(0,0,1,-1,cME1B%P_00_1m1_me,HO%NLEV)
        call Pairs_TMT_JMJ(0,0,1,0,cME1B%P_00_10_me,HO%NLEV)
        call Pairs_TMT_JMJ(0,0,1,1,cME1B%P_00_1p1_me,HO%NLEV)

!       ..... iso-vector pairing: nn,pp,np (T=1; J=0)
        if(.NOT. ALLOCATED(cME1B%P_1m1_00_me))  ALLOCATE(cME1B%P_1m1_00_me(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%P_10_00_me))  ALLOCATE(cME1B%P_10_00_me(1:HO%NLEV*HO%NLEV))
        if(.NOT. ALLOCATED(cME1B%P_1p1_00_me))  ALLOCATE(cME1B%P_1p1_00_me(1:HO%NLEV*HO%NLEV))

!       ...... Only the T=1,J=0 components considered (J can be nonzero ?)
        call Pairs_TMT_JMJ(1,-1,0,0,cME1B%P_1m1_00_me,HO%NLEV)     !pp
        call Pairs_TMT_JMJ(1,0,0,0,cME1B%P_10_00_me,HO%NLEV)       !pn
        call Pairs_TMT_JMJ(1,1,0,0,cME1B%P_1p1_00_me,HO%NLEV)      !nn

        return
        end
!       ..............................................................
        subroutine adjust_cons(NCONS,U0,V0,THoul_work,U1,V1,ALcholes,NLEV)
!       ..............................................................
        USE VAPHFB_PAR
        implicit real*8(a-h,o-z)
        
        !complete matrices
        DIMENSION U0(NLEV,NLEV),V0(NLEV,NLEV)
        DIMENSION U1(NLEV,NLEV),V1(NLEV,NLEV)
        DIMENSION THoul_work(NLEV*NLEV)
        DIMENSION ALcholes(NLEV,NLEV)

        
        real*8, DIMENSION(:), allocatable :: Bmatrix,AMV
        integer, dimension(:), allocatable :: IPIV
        real*8, DIMENSION(:,:), allocatable :: Amatrix
        real*8, DIMENSION(:,:), allocatable :: RO_1,Akapa10_1,Akapa01_1 
  
!        DIMENSION Bmatrix(NCONS)
!        DIMENSION Amatrix(NCONS,NCONS)
!        DIMENSION AMV(NCONS)
!        DIMENSION IPIV(NCONS)

!        DIMENSION RO_1(NLEV,NLEV)
!        DIMENSION Akapa10_1(NLEV,NLEV)
!        DIMENSION Akapa01_1(NLEV,NLEV)

!        dimension ACONS_20(1:HO%NLEV*HO%NLEV*NCONS)
!        dimension AUX1(1:HO%NLEV*HO%NLEV*NCONS*(NCONS+1)/2)

        real*8, DIMENSION(:), allocatable :: ACONS_20 !(NCONS*NLEV*NLEV)
        real*8, DIMENSION(:), allocatable :: AUX1 

        if(.NOT. ALLOCATED(ACONS_20))  ALLOCATE(ACONS_20(1:HO%NLEV*HO%NLEV*NCONS))
        if(.NOT. ALLOCATED(AUX1))  ALLOCATE(AUX1(1:HO%NLEV*HO%NLEV*NCONS*(NCONS+1)/2))


        if(.NOT. ALLOCATED(Bmatrix))  ALLOCATE(Bmatrix(1:NCONS))
        if(.NOT. ALLOCATED(AMV))      ALLOCATE(AMV(1:NCONS))
        if(.NOT. ALLOCATED(IPIV))     ALLOCATE(IPIV(1:NCONS))
        if(.NOT. ALLOCATED(Amatrix))  ALLOCATE(Amatrix(1:NCONS,1:NCONS))
        if(.NOT. ALLOCATED(RO_1))     ALLOCATE(RO_1(1:NLEV,1:NLEV))
        if(.NOT. ALLOCATED(Akapa10_1))     ALLOCATE(Akapa10_1(1:NLEV,1:NLEV))
        if(.NOT. ALLOCATED(Akapa01_1))     ALLOCATE(Akapa01_1(1:NLEV,1:NLEV))

        DO llll=1,Niter_cons  ! 100

        call Lmatrix(THoul_work,ALcholes,NLEV)

        call newUV(U0,V0,THoul_work,ALcholes,U1,V1,NLEV)
         
        call UV2Density(U1,V1,RO_1,Akapa10_1,Akapa01_1,NLEV)

!     computation of AMV (expectation value) from ACONS_ME (ME1B)
!     ................
!     constraints:
!     ACONS_ME(1)                    for  Z
!     ACONS_ME(1+  NLEV**2)          for  N
!     ACONS_ME(1+(2,3,..,7)*NLEV**2) for Q2mu
!     ACONS_ME(1+(8,10,10)*NLEV**2)  for Jx,Jy,Jz
!     ACONS_ME(1+(11,..,16)*NLEV**2) for pairing 
!     ................

        call N_0(V1,AMV(2),AMV(1),NLEV)

!      ............................
!      All the constraints except for N&Z, and pairing: 
!      Q2mu and Jx,Jy,Jz
        do mm=3,NCONS - ipair_const      
           call Q_0(ACONS_ME(1+(mm-1)*NLEV**2),RO_1,Q_0p,Q_0n,NLEV)
           AMV(mm)=Q_0p+Q_0n
        end do        

!      Constraints solely for pairing
        do mm=NCONS-ipair_const+1,NCONS 
           call Qpair_0(ACONS_ME(1+(mm-1)*NLEV**2),Akapa10_1,Akapa01_1,Qpair0,NLEV)
           AMV(mm)=Qpair0
        end do

        
        ia=0
!       ......... All the constraints except for the pairing terms
        do ii=1,NCONS-ipair_const
           indi=(ii-1)*NLEV*NLEV
           Bmatrix(ii)=acons_mv(ii)-AMV(ii)   !  O_i - <O_i>
           if(dabs(Bmatrix(ii)).gt.Input%tolcons) ia = ia+1
            call F_20(U1,V1,ACONS_ME(indi+1),ACONS_20(indi+1),NLEV)
        end do

           acons_mv(3) = qtc20  ! set acons_mv back to original q_20 
!      Constraints solely for pairing
        do ii=NCONS-ipair_const+1,NCONS
           indi=(ii-1)*NLEV*NLEV
           Bmatrix(ii)=acons_mv(ii)-AMV(ii)
           if(dabs(Bmatrix(ii)).gt.Input%tolcons) ia=ia+1
            call G_20(U1,V1,ACONS_ME(indi+1),ACONS_20(indi+1),NLEV)
        end do
        if(ia.eq.0) exit


!!!!!!!!!!!!!!!!!!!!!!!

        lll=0
        do ii=1,NCONS
         indi=(ii-1)*NLEV*NLEV+1
         do jj=ii,NCONS
         indj=(jj-1)*NLEV*NLEV+1
         lll=lll+1
         indl=(lll-1)*NLEV*NLEV+1
         
        call DGEMM ('t','n',NLEV,NLEV,NLEV,one,                 &
     &  ACONS_20(indi),NLEV,ACONS_20(indj),                     &
     &  NLEV,zero,AUX1(indl),NLEV)
         end do
        end do
        
        lll=0
        do ii=1,NCONS
         do jj=ii,NCONS
         lll=lll+1
         indl=(lll-1)*NLEV*NLEV+1

           Amatrix(jj,ii)=zero  
           do mm=1,NLEV
              idig=(mm-1)*NLEV+mm-1
              Amatrix(jj,ii)=Amatrix(jj,ii)+AUX1(indl+idig)
           end do
         
         end do
         
        end do
        
        do ii=1,NCONS
         do jj=ii,NCONS
           Amatrix(ii,jj)=Amatrix(jj,ii)   
         end do
        end do

!       ........... 
!       A * X = B 
!       If INFO=0 (success), then B stores the solution X
!      ............           
        call DGESV(NCONS,1,Amatrix,NCONS,IPIV,Bmatrix,NCONS,INFO)
        
        
        lll=0
        do ii=1,NLEV
         do jj=1,NLEV
         lll=lll+1
         do mm=1,NCONS
            THoul_work(lll)=THoul_work(lll)+                   &
     &      Bmatrix(mm)*ACONS_20(lll+(mm-1)*NLEV*NLEV)
             end do
         end do
        end do
!      ...................            
        END DO  ! llll

        if(ALLOCATED(ACONS_20))  DEALLOCATE(ACONS_20)
        if(ALLOCATED(AUX1))      DEALLOCATE(AUX1)

        if(ALLOCATED(Amatrix))   DEALLOCATE(Amatrix)
        if(ALLOCATED(Bmatrix))   DEALLOCATE(Bmatrix)
        if(ALLOCATED(AMV))       DEALLOCATE(AMV)
        if(ALLOCATED(IPIV))      DEALLOCATE(IPIV)
        if(ALLOCATED(RO_1))      DEALLOCATE(RO_1)
        if(ALLOCATED(Akapa10_1))      DEALLOCATE(Akapa10_1)
        if(ALLOCATED(Akapa01_1))      DEALLOCATE(Akapa01_1)

        end subroutine



!       ...................................1.d0,.............
        subroutine newUV(U0,V0,THoul,ALcholes,U1,V1,NLEV)
!       ....................................................
        USE VAPHFB_PAR
        implicit real*8(a-h,o-z)
        
        !complete matrices
        DIMENSION U0(NLEV,NLEV),V0(NLEV,NLEV)
        DIMENSION U1(NLEV,NLEV),V1(NLEV,NLEV)
        DIMENSION ALcholes(NLEV,NLEV)
        DIMENSION THoul(NLEV,NLEV)
        DIMENSION AUX1(NLEV,NLEV)
        DIMENSION AUX2(NLEV,NLEV)
        DIMENSION AUX3(NLEV,NLEV)
        DIMENSION AUX4(NLEV,NLEV)
        
                
        

        CALL DTRTRI( 'L', 'N', NLEV, ALcholes, NLEV, INFO )
        if(INFO.ne.0) then
         print*,'Error in the inverse'
        end if

        

        call DGEMM ('n','n',NLEV,NLEV,NLEV,1.d0,                   &
     &              V0,NLEV,THoul,NLEV,0.d0,AUX1,NLEV)   !THoul^+*V0^+
                
        call DGEMM ('n','n',NLEV,NLEV,NLEV,1.d0,                   &
     &              U0,NLEV,THoul,NLEV,0.d0,AUX2,NLEV)   !THoul^+*U0^+
         


        do ii=1,NLEV
          do jj=1,NLEV
            AUX3(jj,ii)=U0(jj,ii)+AUX1(jj,ii)	!U0+V0*THoul
            AUX4(jj,ii)=V0(jj,ii)+AUX2(jj,ii)	!V0+U0*THoul

          end do
        end do
       
        call DGEMM ('n','t',NLEV,NLEV,NLEV,1.d0,  &
     &     AUX3,NLEV,ALcholes,NLEV,0.d0,U1,NLEV)  !                    -1t
                                               !U1=(U0+V0*THoul)*L0   

        call DGEMM ('n','t',NLEV,NLEV,NLEV,1.d0,  &
     &     AUX4,NLEV,ALcholes,NLEV,0.d0,V1,NLEV)  !                    -1t

                                               !V1=(UV+U0*THoul)*L0  
        end subroutine

!     ..........................................
