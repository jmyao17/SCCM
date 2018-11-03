       subroutine Prep(lpr)
       USE VAPHFB_Par
       implicit none 

       logical lpr
       real*8  x1,x2

!    .................. PNP
!        PNP%NFOM      = 1
!    .................. AMP
!        AMP%NLEG_ALP  = 1
!        AMP%NLEG_BET  = 1
!        AMP%NLEG_GAM  = 1
        AMP%NNNP      = AMP%NLEG_BET*AMP%NLEG_ALP*AMP%NLEG_GAM
        AMP%JJMax  = JJmax 
!      ...................
!      0: Baranger-Kumar; 
!      1: Standard
!      ...................

       Const%iQBType = 0   
!      ....................      
       H%iden = 0
       H%Aref = 1
       HO%hb0   = hbc/(two*amu)
       HO%hw    = 41.0*Nucl%nucleon(2)**(-third)
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
       endif

!     ............................... SM int.
       if(Input%cIntID  == "GCN" .or. Input%cIntID  == "GCN2850") then
          H%Aref  = 1 
          H%iden  = 0
          HO%emax = 6 !4
          HO%lmax = HO%emax
          HO%twojmax = HO%lmax*2+1
          HO%nmax    = HO%emax/2+1 
          HO%ljmax  = HO%lmax*2 
          HO%nljmax = 4
          HO%NLEV   = 44
          Nucl%ncore(0:2)= (/28,28,56/)  ! (/20,20,40/) (pf)         ! (/8,8,16/) (sd)
       endif

      if(Input%cIntID(1:3)  == "KB3" ) then
          H%Aref  = 42
          H%iden  = 2
          HO%emax = 6 
          HO%lmax = HO%emax 
          HO%twojmax = HO%lmax*2+1  ! twojmax=13
          HO%nmax    = HO%emax/2+1
          HO%ljmax  = HO%lmax*2
          HO%nljmax = 4
          HO%NLEV   = 40
          Nucl%ncore(0:2)= (/20,20,40/)      ! (/8,8,16/) (sd)
       endif

       if(Input%cIntID  == "USD" ) H%Aref=18 

       HO%jmaxp5 = (HO%twojmax+1)/2
!      .............................. mass scale
       if(Input%cIntID .eq. "USD") H%iden = 1

!      ............... input
        Nucl%nucleon(0)= Input%nneut+Nucl%ncore(0) ! neutron number
        Nucl%nucleon(1)= Input%nprot+Nucl%ncore(1) ! neutron number
        Nucl%nucleon(2)= Nucl%nucleon(0)+Nucl%nucleon(1) ! total

        Nucl1%nucleon(0) = Nucl%nucleon(0)-2
        Nucl1%nucleon(1) = Nucl%nucleon(1)+2
        Nucl1%nucleon(2) = Nucl%nucleon(2)
!     ................ basis
        HO%b_osc = sqrt(hbc*two*HO%hb0/HO%hw)

        if(lpr) then
         write(*,*) ' Oscillator Length of HO: ',HO%b_osc
         write(*,'(a,3i6)') ' # of Nucleons in Core:',Nucl%ncore
         write(*,'(a,3i6)') '   Total # of Nucleons:',Nucl%nucleon
        endif
      call nucleus(1,Nucl%nucleon(1),Nucl%nucnam)
      call nucleus(1,Nucl1%nucleon(1),Nucl1%nucnam)

!     conefficient between beta and Q value
!      beta_L = 4*pi/(3AR^L) * Q_L
!     .......................................
       HO%hw = 41.20*Nucl%nucleon(2)**(-1./3.)

       if(Const%iQBType.eq.0) then
          Const%Q2BN = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BP = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
          Const%Q2BA = 0.4d0/HO%hw  ! 11.33661688d0          ! KB
          print*, ' Const%Q2BA=', Const%Q2BA
       else
          Const%Q2BA=(4.*pi)/(3.*(r0**2)*(Nucl%nucleon(2)**(5./3.))) !factor for quadrupole moments
          Const%Q2BP=(4.*pi)/(3.*(r0**2)*Nucl%nucleon(1)*(Nucl%nucleon(2)**(2./3.)))
          Const%Q2BN=(4.*pi)/(3.*(r0**2)*Nucl%nucleon(0)*(Nucl%nucleon(2)**(2./3.)))
       endif

       x1= Nucl%nucleon(2) ! mass number A = xp1+xn1
       x2= H%Aref         ! xp2+xn2    ! 18

!    ....................... parameter for two-body states
       TPB%bJMax= 2*HO%lmax + 1

       return
       end

