       subroutine Prep(lpr)
       USE VAPHFB_Par
       implicit none 
       logical lpr  
       real*8  x1,x2

!    ....... truncation model space with a cutoff
!    .................. PNP
!        PNP%NFOM      = 1
!    .................. AMP
!        AMP%NLEG_ALP  = 1
!        AMP%NLEG_BET  = 1
!        AMP%NLEG_GAM  = 1
        AMP%NNNP      = AMP%NLEG_BET*AMP%NLEG_ALP*AMP%NLEG_GAM
        AMP%JJMax  = 6

          Input%en=0.0
          Input%ep=1.0
!      ...................
!      0: Baranger-Kumar; 
!      1: Standard
!      ...................

       Const%iQBType = 0   
!      ....................      
       H%iden = 0
       H%Aref = 1
       HO%hb0   = hbc/(two*amu)
!     ............................... Chiral int.
          HO%tmax = 1
       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax04' .or. Input%cValID == 'eMax04' )) then
!          print *, ' Chiral Interaction Is Adopted '
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
       endif

       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax06' .or. Input%cValID == 'eMax06' )) then
          if(lpr) print *, ' Chiral Interaction Is Adopted '
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
       endif
       if(Input%IntType .eq. 1 .and. (Input%cValID == 'emax08' .or. Input%cValID == 'eMax08' )) then
!          print *, ' Chiral Interaction Is Adopted '
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

!     ............................... SM int.
       if(Input%IntType .eq. 0 .and. Input%cIntID(1:3)  == "GCN" ) then
          Input%en=0.5
          Input%ep=1.5
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

      if(Input%IntType .eq. 0 .and. Input%cIntID(1:3)  == "KB3" ) then
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
       endif

       HO%jmaxp5 =  (HO%twojmax+1)/2 ! emax=4, twojmax = 9, HO%jmaxp5=5
       HO%tnljmax = HO%nljmax*2 

!      ............... input
        Nucl%nucleon(0)= Input%nneut+Nucl%ncore(0) ! neutron number
        Nucl%nucleon(1)= Input%nprot+Nucl%ncore(1) ! neutron number
        Nucl%nucleon(2)= Nucl%nucleon(0)+Nucl%nucleon(1) ! total

        if(lpr) then
        print *,' -> Information about the Nucleus: '
        write(*,'(a,i3)') '     Proton(v): ',Input%nprot
        write(*,'(a,i3)') '     Proton(c): ',Nucl%ncore(1)
        write(*,'(a,i3)') '     Proton(t): ',Nucl%nucleon(1)
        write(*,'(a,i3)') '     Neutron(v):',Input%nneut
        write(*,'(a,i3)') '     Neutron(c):',Nucl%ncore(0)
        write(*,'(a,i3)') '     Neutron(t):',Nucl%nucleon(0)
        endif
!        HO%b_osc = 1.01 * Nucl%nucleon(2)**(1.d0/6.d0)  ! ring & schuck's book (2.13)
       
      call nucleus(1,Nucl%nucleon(1),Nucl%nucnam)
!     conefficient between beta and Q value
!      beta_L = 4*pi/(3AR^L) * Q_L
!     .......................................
       HO%hw = 41.20*Nucl%nucleon(2)**(-1./3.)
       if(Input%IntType .eq. 1 ) then
           HO%hw    = Input%ihwHO ! frequency of the HO basis for the chiral interaciton 
       else
           HO%hw = 41.20*Nucl%nucleon(2)**(-1./3.)
           Input%ihwHO = Int(HO%hw)
       endif
       HO%b_osc = sqrt(hbc*two*HO%hb0/HO%hw)
       if(lpr) then
        print *,   ' -> Information about the HO basis: '
        print *,   '    Frequency of HO: ',HO%hw
        write(*,*) '    Oscillator Length of HO: ',HO%b_osc
       endif
       if(Const%iQBType.eq.0) then
          Const%Q2BN = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BP = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BA = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
!          print*, ' Const%Q2BA=', Const%Q2BA
       else
          Const%Q2BA=(4.*pi)/(3.*(r0**2)*(Nucl%nucleon(2)**(5./3.))) !factor for quadrupole moments
          Const%Q2BP=(4.*pi)/(3.*(r0**2)*Nucl%nucleon(1)*(Nucl%nucleon(2)**(2./3.)))
          Const%Q2BN=(4.*pi)/(3.*(r0**2)*Nucl%nucleon(0)*(Nucl%nucleon(2)**(2./3.)))
       endif

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
        else if (H%iden.eq.3) then
         H%ddd_kin=(x2/x1)**(0.3d0)
         H%ddd_tbme=(x2/x1)**(0.3d0)
        end if
!    ....................... parameter for two-body states
       TPB%bJMax= 2*HO%lmax + 1
       return
       end

