      subroutine Filename4collwf(NOQ)
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
      character*6 cname0

      cname0 = '../../'
      name_NOQ1  =mod(NOQ/10,10)+48
      name_NOQ2  =mod(NOQ,10)+48
      ihw1= mod(Input%ihwHO/10,10)+48
      ihw2= mod(Input%ihwHO,10)+48

      File%FF=cname0//Nucl%nucnam//'/F4GS.'         &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //'NOQ'//char(name_NOQ1)//char(name_NOQ2)//'.dat'


      end
!.......................................................................
      subroutine Filename4Kernels(icase,iq1,iq2)
!.......................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
!      implicit none
      INTEGER nuc1,nuc2,iq1,iq2,ihw1,ihw2
      character*1 sign1,sign2
      character*6 cname0
!---------------------------
              if(icase.eq.0) name   = 1+48
              if(icase.eq.1) name   = 3+48
              betac1 = Const%beta2t_mesh(iq1)
              betac2 = Const%beta2t_mesh(iq2)
              gammac1= Const%gamma2t_mesh(iq1)
              gammac2= Const%gamma2t_mesh(iq2)
              r2c1   = Const%P00_mesh(iq1) 
              r2c2   = Const%P00_mesh(iq2)
              cf1    = Const%hw_mesh(iq1)
              cf2    = Const%hw_mesh(iq2)
              cname0 = '../../'
              if(betac1.ge.0.d0) sign1='+'
              if(betac1.lt.0.d0) sign1='-'
              if(betac2.ge.0.d0) sign2='+'
              if(betac2.lt.0.d0) sign2='-'
              ab2c1  = abs(betac1)
              name11 = ab2c1+48
              name21 = mod(ab2c1*10,10.d0)  +48
              name31 = mod(ab2c1*100,10.d0) +48
              name41 = mod(gammac1/10,10.d0)+48
              name51 = mod(gammac1,10.d0)   +48

              name61 = mod(Int(r2c1*1.001),10)+48
              name71 = mod(Int(r2c1*1.001)/10,10) +48
              namer1 = mod(r2c1*100.0001,10.d0) + 48

              ncf11 = cf1+48
              ncf12 = mod(cf1*10,10.d0)  +48
              ncf13 = mod(cf1*100,10.d0) +48
!              name61 = r2c1+48
!              name71 = mod(r2c1*10,10.d0) +48
!              namer1 = mod(r2c1*100.0001,10.d0) +48
!---------------------------
              ab2c2  = abs(betac2)
              name12 = ab2c2+48
              name22 = mod(ab2c2*10,10.d0)  +48
              name32 = mod(ab2c2*100,10.d0) +48
              name42 = mod(gammac2/10,10.d0)+48
              name52 = mod(gammac2,10.d0)   +48

              name62 = mod(Int(r2c2*1.001),10)+48
              name72 = mod(Int(r2c2*1.001)/10,10) +48
              namer2 = mod(r2c2*100.0001,10.d0) + 48

              ncf21 = cf2+48
              ncf22 = mod(cf2*10,10.d0)  +48
              ncf23 = mod(cf2*100,10.d0) +48
!              namep1 = mod(Int(p00*1.001)/10,10)  +48
!              namep2 = mod(Int(p00*1.001),10)  +48
!---------------------------  
              jphi    = PNP%MPhi !PNP%NFOM 
              name81  = mod(jphi/10,10) + 48
              name82  = mod(jphi,10) + 48
              name91  = mod(AMP%NLEG_BET/10,10) + 48
              name92  = mod(AMP%NLEG_BET,10) + 48

              nalp1   = mod(AMP%NLEG_ALP/10,10) + 48
              nalp2   = mod(AMP%NLEG_ALP,10) + 48
              ngam1   = mod(AMP%NLEG_GAM/10,10) + 48
              ngam2   = mod(AMP%NLEG_GAM,10) + 48

              name_emax1 =mod(HO%emax/10,10) + 48               
              name_emax2 =mod(HO%emax,10) + 48               

              ihw1= mod(Input%ihwHO/10,10)+48
              ihw2= mod(Input%ihwHO,10)+48


              nuc1 = mod(Nucl%nucleon(2)/10,10) +48
              nuc2 = mod(Nucl%nucleon(2),10)    +48
!              if(mphi(1).lt.10)
!     &         name6  = '0'//char(jphi)
!---------------------------  
      File%elem =cname0//Nucl%nucnam//'/kern.'//char(name)//'D'         &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(nalp1)//char(nalp2)       &
     &          //'.'//char(name91)//char(name92)      &
     &          //'.'//char(ngam1)//char(ngam2)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&  
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.elem'
      File%Rho1B=cname0//Nucl%nucnam//'/Rho1B.'//Nucl%nucnam//char(nuc1)//char(nuc2) &
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(nalp1)//char(nalp2)       &
     &          //'.'//char(name91)//char(name92)      &
     &          //'.'//char(ngam1)//char(ngam2)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'

      File%TD1B=cname0//Nucl%nucnam//'/TD1B.'//Nucl%nucnam//char(nuc1)//char(nuc2) &
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(nalp1)//char(nalp2)      &
     &          //'.'//char(name91)//char(name92)       &
     &          //'.'//char(ngam1)//char(ngam2)//'.'       &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'

      File%Rho2B=cname0//Nucl%nucnam//'/Rho2B.'//Nucl%nucnam//char(nuc1)//char(nuc2) &
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(nalp1)//char(nalp2)       &
     &          //'.'//char(name91)//char(name92)      &
     &          //'.'//char(ngam1)//char(ngam2)//'.'       &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'


      File%Rho3B=cname0//Nucl%nucnam//'/Rho3B.'//Nucl%nucnam//char(nuc1)//char(nuc2) &
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(nalp1)//char(nalp2)       &
     &          //'.'//char(name91)//char(name92)       &
     &          //'.'//char(ngam1)//char(ngam2)//'.'       &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(ihw1)//char(ihw2)//'_' &
     &           //trim(Input%cFlow)//'_'                      &
     &          //trim(Input%vs4me3b)//'_' &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'

!--------------------- 
       ! print *, File%Rho1B,File%Rho2B
      return
      end
!________________________________________________________________________
      subroutine Filename4wfs(betat,gammat,p00,cf)
!.......................................................................
      USE VAPHFB_PAR
      implicit none
      character*1 sign1,sign2
      character*6 cname0
      real*8      betat,gammat,p00,ab2c1,cf
      integer     name_emax1,name_emax2,nuc1,nuc2,name1,name2,name3,&
     &            name4,name5,namep1,namep2,namep3,hw1,hw2,cf1,cf2,cf3
!---------------------------
              cname0 = '../../'
              if(betat.ge.0.d0) sign1='+'
              if(betat.lt.0.d0) sign1='-'
              ab2c1  = abs(betat)
              name1 = ab2c1+48
              name2 = mod(ab2c1*10,10.d0)  +48
              name3 = mod(ab2c1*100,10.d0) +48
              name4 = mod(gammat/10,10.d0)+48
              name5 = mod(gammat,10.d0)   +48
              nuc1 = mod(Nucl%nucleon(2)/10,10) +48
              nuc2 = mod(Nucl%nucleon(2),10)    +48

!              namep1 = p00+48
!              namep2 = mod(p00*10,10.d0)  +48
!              namep3 = mod(p00*100,10.d0) +48

              namep1 = mod(Int(p00*1.001)/10,10)  +48
              namep2 = mod(Int(p00*1.001),10)  +48

              name_emax1 =mod(HO%emax/10,10) + 48
              name_emax2 =mod(HO%emax,10) + 48

              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48

              cf1 = cf+48
              cf2 = mod(cf*10,10.d0)  +48
              cf3 = mod(cf*100,10.d0) +48
!---------------------------  
       if(Input%IsHFB.eq.0) then

      File%wf =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_HF_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)                 &
     &          //'cf'//char(cf1)//char(cf2)//char(cf3)// '.wf' !  

      else if(Input%IsHFB.eq.2) then
      File%wf =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_VAP_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)                 &
     &          //'cf'//char(cf1)//char(cf2)//char(cf3)// '.wf' ! 

      else if(Input%IsHFB.eq.1) then
      File%wf =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_HFB_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)                 &
     &          //'cf'//char(cf1)//char(cf2)//char(cf3)// '.wf' ! 
      endif

!---------------------------  
      if(Input%IsHFB.eq.2) then
      File%out =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_VAP_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      

      else if(Input%IsHFB.eq.1) then
      File%out =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_HFB_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &        //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      

      else if(Input%IsHFB.eq.0) then
      File%out =cname0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_HF_'//trim(Input%cIntID)       &
!     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &        //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      
      endif
     
      return
      end





