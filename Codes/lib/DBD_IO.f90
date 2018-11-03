!.......................................................................
      subroutine Filename4cwfs()
!.......................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
!      implicit none
      INTEGER iq1,iq2,hw1,hw2
      character*1 sign1,sign2
      character*6 name0
!---------------------------
              if(icase.eq.0) name   = 1+48
              if(icase.eq.1) name   = 3+48
             
              name0 = '../../'
              name_emax1 =mod(HO%emax/10,10) + 48               
              name_emax2 =mod(HO%emax,10) + 48               
              name_NOQi1  =mod(Input%NGCM2/10,10)+48
              name_NOQi2  =mod(Input%NGCM2,10)+48
              name_NOQf1  =mod(Input%NGCM1/10,10)+48
              name_NOQf2  =mod(Input%NGCM1,10)+48
              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48
!--------------------- 
      File%FFi=name0//Nucl%nucnam//'/F4GS.'         &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &        //trim(Input%cFlow)//'_'                       &
     &          //'NOQ'//char(name_NOQi1)//char(name_NOQi2)//'.dat'

      File%FFf=name0//Nucl1%nucnam//'/F4GS.'         &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &        //trim(Input%cFlow)//'_'                       &
     &          //'NOQ'//char(name_NOQf1)//char(name_NOQf2)//'.dat'
!--------------------- 
      return
      end
!______________________________________________________________________________
!.......................................................................
      subroutine Filename4Kernels(icase,iq1,iq2)
!.......................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
!      implicit none
      INTEGER hw1,hw2,iq1,iq2
      character*1 sign1,sign2
      character*6 name0
!---------------------------
              if(icase.eq.0) name   = 1+48
              if(icase.eq.1) name   = 3+48
              betac1 = Const1%beta2t_mesh(iq1)
              betac2 = Const%beta2t_mesh(iq2)
              gammac1= Const1%gamma2t_mesh(iq1)
              gammac2= Const%gamma2t_mesh(iq2)
              r2c1   = Const1%P00_mesh(iq1) 
              r2c2   = Const%P00_mesh(iq2)
             
              name0 = '../../'
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

!              name62 = r2c2+48
!              name72 = mod(r2c2*10,10.d0) +48
!              namer2 = mod(r2c2*100.0001,10.d0) +48
!---------------------------  
              jphi    = PNP%NFOM 
              name81  = mod(jphi/10,10) + 48
              name82  = mod(jphi,10) + 48
              name91  = mod(AMP%NLEG_BET/10,10) + 48
              name92  = mod(AMP%NLEG_BET,10) + 48

              name_emax1 =mod(HO%emax/10,10) + 48               
              name_emax2 =mod(HO%emax,10) + 48               


              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48
!--------------------- 
      File%TD2B =name0//Nucl%nucnam//'/TD2B.'//char(name)//'D'     &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //trim(Input%cIntID)//'_'   &
     &           //trim(Input%cFlow)//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)//    &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)//'.dens'

!--------------------- 
      return
      end

!.......................................................................
      subroutine Filename4Results()
!.......................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
!      implicit none
      INTEGER hw1,hw2,iq1,iq2
      character*1 sign1,sign2
      character*6 name0
!---------------------------
              jphi    = PNP%NFOM 
              name81  = mod(jphi/10,10) + 48
              name82  = mod(jphi,10) + 48
              name91  = mod(AMP%NLEG_BET/10,10) + 48
              name92  = mod(AMP%NLEG_BET,10) + 48

              name_emax1 =mod(HO%emax/10,10) + 48               
              name_emax2 =mod(HO%emax,10) + 48               

              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48
!--------------------- 
      File%out='M0nu_'//Nucl%nucnam//'2'//Nucl1%nucnam &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &          //trim(Input%cIntID)//'_'       &
     &          //trim(Input%cValID)   &
     &          //'.out'
      File%tdens='M0nu_'//Nucl%nucnam//'2'//Nucl1%nucnam &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
!     &          //'eMax'//char(name_emax1)//char(name_emax2) &
     &          //trim(Input%cIntID)//'_'       &
     &          //trim(Input%cValID)   &
     &          //'.tdens'

      File%imsrg='M0nu_'//Nucl%nucnam//'2'//Nucl1%nucnam &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
!     &          //'eMax'//char(name_emax1)//char(name_emax2) &
     &          //trim(Input%cIntID)//'_'       &
     &          //trim(Input%cValID)   &
     &          //'.tdens4imsrg'
!--------------------- 
      return
      end
