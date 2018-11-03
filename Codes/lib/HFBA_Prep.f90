
       subroutine Prep(lpr)

       USE VAPHFB_Par
       implicit none 
       logical lpr
       real*8 chi,R00


       real*8  x1,x2

      ! parameters for optimizers
      !Opt%mp   =0.25 
      !Opt%bet1 =0.33 !0.9 
      !Opt%bet2 =0.d0 !0.99 < 1 
      Opt%eps=1.d-8 
      Opt%gam1 = 0.75 !0.90  ! 0.75  bet1/(1-bet1)
      Opt%gam2 = 0.d0  ! bet2/(1-bet2)

      HO%tp(1) = ' +'
      HO%tp(2) = ' -'
      HO%tis(1) = 'n'
      HO%tis(2) = 'p'
      HO%tl(0)  = 's'
      HO%tl(1)  = 'p'
      HO%tl(2)  = 'd'
      HO%tl(3)  = 'f'
      HO%tl(4)  = 'g'
      HO%tl(5)  = 'h'
      HO%tl(6)  = 'i'
      HO%tl(7)  = 'j'
      HO%tl(8)  = 'k'
      HO%tl(9)  = 'l'

!    .................. PNP
!        PNP%NFOM      = 1
!    .................. AMP
!        AMP%NLEG_ALP  = 1
!        AMP%NLEG_BET  = 8
!        AMP%NLEG_GAM  = 1
!      ...................
!      0: Baranger-Kumar; 
!      1: Standard
!      ...................
       H%iabcd_max   = 0
       Const%iQBType = 0   
       chi   = 0.60
!      ....................      
       H%iden = 0
       H%Aref = 1
       HO%hb0   = hbc/(two*amu)
!       HO%hw    = 41.0*Nucl%nucleon(2)**(-third)
!     ............................... Chiral int.
          HO%tmax = 1
       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax04' .or. Input%cValID == 'eMax04' )) then
!          print *, ' Chiral Interaction Is Adopted '
!          call Model_Space_Generator(4,HO%nljmax)
          H%Aref  = 1
          H%iden  = 0
          HO%emax = 4 !4
          HO%lmax = HO%emax 
          HO%twojMax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2 
          HO%nljmax = 15 
          HO%NLEV   = 140 
          HO%NLEV2  = HO%NLEV**2
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Const%iQBType = 1  !   
          HO%hw    = Input%ihwHO ! frequency of the HO basis for the chiral interaciton 
       endif
       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax06' .or. Input%cValID == 'eMax06' )) then
!          print *, ' Chiral Interaction Is Adopted '
          H%Aref  = 1
          H%iden  = 0
          HO%emax = 6 !4
          HO%lmax = HO%emax 
          HO%twojMax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2 
          HO%nljmax = 28 
          HO%NLEV   = 336 
          HO%NLEV2  = HO%NLEV**2
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Const%iQBType = 1  !   
          HO%hw    = Input%ihwHO ! frequency of the HO basis for the chiral interaciton 
       endif

       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax08' .or. Input%cValID == 'eMax08' )) then
!          print *, ' Chiral Interaction Is Adopted '
!          call Model_Space_Generator(8,HO%nljmax)
          H%Aref  = 1
          H%iden  = 0
          HO%emax = 8 !4
          HO%lmax = HO%emax
          HO%twojMax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2
          HO%nljmax = 45 
          HO%NLEV   = 660 
          HO%NLEV2  = HO%NLEV**2
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Const%iQBType = 1  !   
          HO%hw    = Input%ihwHO ! frequency of the HO basis for the chiral interaciton 
       endif

       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax10' .or. Input%cValID == 'eMax10' )) then
!          print *, ' Chiral Interaction Is Adopted '
!          call Model_Space_Generator(10,HO%nljmax)
          H%Aref  = 1
          H%iden  = 0
          HO%emax = 10
          HO%lmax = HO%emax
          HO%twojMax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2
          HO%nljmax = 66
          HO%NLEV   = 1144 
          HO%NLEV2  = HO%NLEV**2
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Nucl%ncore(0:2)= (/0,0,0/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          Const%iQBType = 1  !   
          HO%hw    = Input%ihwHO ! frequency of the HO basis for the chiral interaciton 
       endif

!     ............................... SM int.
       if(Input%cIntID(1:3)  == "GCN" ) then
          H%Aref  = 1 
          H%iden  = 0
          HO%emax = 6 !4
          HO%lmax = HO%emax !4
          HO%twojmax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1 
          HO%ljmax  = HO%lmax*2 
          HO%nljmax = 4
          HO%NLEV   = 44
          Nucl%ncore(0:2)= (/28,28,56/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
          chi   = 0.28
       endif

       if(Input%IntType .eq. 0 .and. Input%cIntID(1:3)  == "USD" ) then 
          H%Aref  = 18 
          H%iden  = 1 ! 1
          HO%emax = 6 !4
          HO%lmax = HO%emax !4
          HO%twojmax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2
          HO%nljmax = 3
          HO%NLEV   = 24 
          Nucl%ncore(0:2)= (/8,8,16/)      ! (/8,8,16/) (sd)
       endif

       if(Input%cIntID(1:3)  == "KB3" ) then
          H%Aref  = 42 
          H%iden  = 2
          HO%emax = 6 !4
          HO%lmax = HO%emax !4
          HO%twojmax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2
          HO%nljmax = 4
          HO%NLEV   = 40
          Nucl%ncore(0:2)= (/20,20,40/)      ! (/8,8,16/) (sd)
          chi   = 0.60
       endif



       HO%jmaxp5 =  (HO%twojmax+1)/2 ! emax=4, twojmax = 9, HO%jmaxp5=5
       HO%tnljmax = HO%nljmax*2

!      ............... input
        Nucl%nucleon(0)= Input%nneut+Nucl%ncore(0) ! neutron number
        Nucl%nucleon(1)= Input%nprot+Nucl%ncore(1) ! proton number
        Nucl%nucleon(2)= Nucl%nucleon(0)+Nucl%nucleon(1) ! total

!       HO%b_osc = 1.01 * Nucl%nucleon(2)**(1.d0/6.d0)  ! ring & schuck's book (2.13)
     
        if(lpr) then 
        print *,' -> Information about the Nucleus: '
        write(*,'(a,i3)') '     Proton(v): ',Input%nprot
        write(*,'(a,i3)') '     Proton(c): ',Nucl%ncore(1)
        write(*,'(a,i3)') '     Proton(t): ',Nucl%nucleon(1)
        write(*,'(a,i3)') '     Neutron(v):',Input%nneut
        write(*,'(a,i3)') '     Neutron(c):',Nucl%ncore(0)
        write(*,'(a,i3)') '     Neutron(t):',Nucl%nucleon(0)
        endif
      call nucleus(1,Nucl%nucleon(1),Nucl%nucnam)

!     conefficient between beta and Q value
!      beta_L = 4*pi/(3AR^L) * Q_L
!     .......................................
!     ................ basis
       if(Input%IntType .eq. 0) then
          HO%hw = 41.20*Nucl%nucleon(2)**(-1./3.)
          Input%ihwHO = Int(HO%hw)
       endif
       HO%b_osc = sqrt(hbc*two*HO%hb0/HO%hw)
       if(lpr) then 
         write(*,*) '   Frequency of HO: ',HO%hw
         write(*,*) '   Oscillator Length of HO: ',HO%b_osc
       endif
       R00   = r0*Nucl%nucleon(2)**(third)  ! R_0
       if(Const%iQBType.eq.0) then
          Const%Q2BN = chi/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BP = chi/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BA = chi/HO%hw  ! 11.33661688d0          ! KB
          print*, ' Const%Q2BA=', Const%Q2BA
       else
          Const%Q2BA=(4.*pi)/(3.d0*R00**2*Nucl%nucleon(2)) ! factor for quadrupole moments
          Const%Q2BP=(4.*pi)/(3.d0*R00**2*Nucl%nucleon(1)) ! factor for quadrupole moments
          Const%Q2BN=(4.*pi)/(3.d0*R00**2*Nucl%nucleon(0)) ! factor for quadrupole moments
       endif

!       if the constraint is placed on the beta instead of Q value
        if(Const%iQB.eq.1) then
           beta2=ACONS_MV_aux(3)
           gamma2=ACONS_MV_aux(4)*pi/180.
           ACONS_MV_aux(3)=beta2*dcos(gamma2)/Const%Q2BA
           ACONS_MV_aux(4)=beta2*dsin(gamma2)/(Const%Q2BA*sqrt(2.d0))  ! Q22

           beta2_IV  = ACONS_MV_aux(6)
           gamma2_IV = ACONS_MV_aux(7)*pi/180.
           ACONS_MV_aux(6)=beta2_IV*dcos(gamma2_IV)/Const%Q2BA
           ACONS_MV_aux(7)=beta2_IV*dsin(gamma2_IV)/(Const%Q2BA*sqrt(2.d0))
        end if

!     ........... scale for the SM 2BME
       x1= Nucl%nucleon(2) ! mass number A = xp1+xn1
       x2= H%Aref         ! xp2+xn2    ! 18

        if (H%iden.eq.0) then
         H%ddd_kin =1.d0
         H%ddd_tbme=1.d0

        else if (H%iden.eq.1) then
         H%ddd_kin=1.d0
         H%ddd_tbme=(x2/x1)**(0.3d0)

        else if (H%iden.eq.2) then
         H%ddd_kin=1.d0
         H%ddd_tbme=(x2/x1)**(third)
         print *, 'H%ddd_tbme=',H%ddd_tbme
        else if (H%iden.eq.3) then
         H%ddd_kin=(x2/x1)**(0.3d0)
         H%ddd_tbme=(x2/x1)**(0.3d0)
        end if

       return
       end

!     ..........................................
      SUBROUTINE ZCONST
!     ..........................................
      USE VAPHFB_PAR
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      COMMON IOUT,IREAD,IPRI
      COMMON /CHAT/ ZHAT(0:200)
      COMMON /CLOG/ ZLOG(100),ZGAMM2(100)
      COMMON /CONST/ ZAC,ZACM,ZPI,ZEPS
      ZAC=SQRT(2.0D0)
      ZACM=1.0D0/ZAC
      ZPI=4.0D0*DATAN(1.0D0)
      ZEPS=1.0D-12
      ZLOG(1)=0.0D0
      ZLOG(2)=0.0D0
      WN=1.0D0
      DO 1 I=3,100
      WN=WN+1.0D0
1     ZLOG(I)=ZLOG(I-1)+DLOG(WN)
      zgamm2(1)=dsqrt(zpi)
      zgamm2(2)=1.0d0
      do jj=3,100
         zgamm2(jj)=(jj-2)*zgamm2(jj-2)/2.0d0
      end do
      DO 3 S=0,200
3     ZHAT(S)=DSQRT(DBLE(S+1))
      RETURN
      END
