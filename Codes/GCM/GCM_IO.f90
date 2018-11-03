!    ............................................
      subroutine table_general(jproj,iexst)
!    ............................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
      real*8 qtot(GCM%NOQ)

  202 format (4i3,4f10.3,2f12.4,3f9.4)

!     ................................. initialization
      !jqk(j,iq,k) = k+j+1+(2*j+1)*(iq-1) 
      maxmp = GCM%qkmax(jproj)
      !print *, 'J=',jproj, 'qkmax=',maxmp
      amas  = Nucl%nucleon(2) 
      r00 = r0*amas**third
      fac = dsqrt(16*pi/5)
!     ................................ print diagonal (in q space)
!     matrix elements
      do iqk1 = 1,maxmp
           !iqk1 = jqk(jproj,iq1,k1) 
           iq1 = GCM%iq(iqk1,jproj)
           k1 = GCM%ik(iqk1,jproj)
      do iqk2 = 1,maxmp
           !iqk2 = jqk(jproj,iq1,k2)
           iq2 = GCM%iq(iqk2,jproj)
           k2  = GCM%ik(iqk2,jproj)
           if(iq1.ne.iq2 .or. k1.ne.k2) cycle    
           ! only print out diagonal elements in q space
           ee  = GCM%hjkkqq(jproj,iqk1,iqk2)
           gg  = GCM%njkkqq(jproj,iqk1,iqk2)
           xnn = dreal(Kernel%nn(iq1,iq1,jproj,k1,k2))
           xj2 = dreal(Kernel%J2(iq1,iq1,jproj,k1,k2))
           xpp = dreal(Kernel%zz(iq1,iq1,jproj,k1,k2))
           qtot(iq1) = (3*amas*r00**2)/(4*pi)*Const%beta2t_mesh(iq1)*fac
           if(abs(gg).gt.CHOP)  &
     &     write(*,202),jproj,k1,k2,iq1,Const%beta2t_mesh(iq1),  &
     &     Const%gamma2t_mesh(iq1),qtot(iq1),Const%P00_mesh(iq1),ee/gg,gg,xnn,xpp,xj2
      enddo ! iq2
      enddo ! iq1
!     ................................ print matrix of existing matrix
!     elements
      do iqk1 = 1,maxmp
      do iqk2 = 1,maxmp
           gg  = GCM%njkkqq(jproj,iqk1,iqk2)
          if(abs(gg).gt.CHOP) then
            GCM%Iexst(iqk1,iqk2) = log(abs(gg))
          endif
      enddo ! iq1
      enddo ! iq2

      return
      end

!______________________________________________________________________________
      subroutine readkernel_general(iq1,iq2,jprmi,jprma,jdf)

!..............................................................................
!     read matrix elements for Hill-Wheeler-Griffin solver to file lou
!     .
!..............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-z)

  401 format (10i5)
  402 format (8e15.8)

!    ......................................................... 
      icheck = 0
      do jpr = jprmi,jprma, jdf
        if(jpr.eq.1) cycle
!       .................................... write matrix elements to
        read (lou,401,end=403) jproj
        if(jproj.ne.jpr) stop 'Error: in reading kernels !'
        read (lou,402) bet1,gam1,bet2,gam2
!    .........................................................
        if(iv(jproj).gt.zero) then         ! for even spin
           jeff    = jproj/2+1
           nmaxdi  = jeff*maxmp
           k1max  =  jproj
        endif
        if(iv(jproj).lt.zero) then       ! for odd spin
           jeff    = (jproj-1)/2
           nmaxdi  = jeff*maxmp
           k1max  =  abs(jproj-1)
        endif
!    .........................................................

         do k1  =  -k1max, k1max, 1
         do k2  =  -k1max, k1max, 1
               read (lou,401) k1r,k2r
               read (lou,402) Kernel%njkk(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%hjkk(iq1,iq2,jproj,k1,k2)

               read (lou,402) Kernel%nn(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%zz(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%J2(iq1,iq2,jproj,k1,k2)

               if(jproj.eq.5) &
              !print*,k1r,k2r,Kernel%nn(iq1,iq2,jproj,k1,k2),Kernel%zz(iq1,iq2,jproj,k1,k2)
              ! read (lou,402) Kernel%qp(iq1,iq2,jproj,k1,k2), Kernel%qcp(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,0),Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,1),Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
               read (lou,402) Kernel%q0p(iq1,iq2,jproj,k1,k2),Kernel%q0pc(iq1,iq2,jproj,k1,k2)
          enddo ! k2
              ! read(lou,402) Kernel%qpred(k1,k1max+2,jproj),Kernel%qcpred(k1,k1max+2,jproj)

               !if(jproj.eq.2) print *,'Kernel%qpred
               !<2||Q2||2>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.4) print *,'Kernel%qpred
               !<4||Q2||4>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.6) print *,'Kernel%qpred
               !<6||Q2||6>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.0) print *,'Kernel%qpred <2||Q2||0>=',
               !Kernel%qpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.0) print
               !*,'Kernel%qcpred<2||Q2||0>=',Kernel%qcpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.2) print *,'Kernel%qpred <4||Q2||2>=',
               !Kernel%qpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.4) print *,'Kernel%qpred <6||Q2||4>=',
               !Kernel%qpred(iq1,iq2,jproj,k1,k2)
!     ....................................... qcp, qcpred were set to
!     zero for iq1=iq2 cases in PKC code
!          q0p    = <J(f),K1,q1 ||e*r^2|| J(i),K2,q2>
!          qp     = <J(f),K1,q1 ||e*Q_2|| J(i),K2,q2>
!          qcp    = <J,K1,q2   ||e*Q_2|| J,K2,q1>
!          qpred  = <J+2,K1,q1 ||e*Q_2|| J,K2,q2>
!          qcpred = <J+2,K1,q2 ||e*Q_2|| J,K2,q1>
!     ........................................ 
         enddo ! k1                           
      icheck = icheck + 1
  403 continue
  99  enddo ! jproj

      if(icheck.eq.0)  write(12,*) elem, 'does not exist ...'

      end 


!    ............................................
      subroutine table(k1,k2,k1m,k2m,jproj,iexst)
!    ............................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
      real*8 qtot(GCM%NOQ)

  202 format (4i3,4f10.3,2f12.4,3f9.4)
      !write(*,201)
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
           if(abs(gg).gt.CHOP)  &
     &     write(*,202),jproj,k1,k2,iq1,Const%beta2t_mesh(iq1),  &
     &               Const%gamma2t_mesh(iq1),qtot(iq1),Const%P00_mesh(iq1),ee/gg,gg,xnn,xpp,xj2
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
!  402 format (8e20.8)

!    ......................................................... 
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
              k1 = 0
              k2 = 0

               read (lou,401) k1r,k2r
               read (lou,402) Kernel%njkk(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%hjkk(iq1,iq2,jproj,k1,k2)

               ! the following line only for temp use, should be deleted  
               !Kernel%hjkk(iq1,iq2,jproj,k1,k2) &
               != Kernel%hjkk(iq1,iq2,jproj,k1,k2)*Kernel%njkk(iq1,iq2,jproj,k1,k2)

               read (lou,402) Kernel%nn(iq1,iq2,jproj,k1,k2), &
     &                        Kernel%zz(iq1,iq2,jproj,k1,k2)

               if(abs(Nucl%nucleon(1)-Kernel%zz(iq1,iq2,jproj,k1,k2)).gt.0.1) &
     &         write(*,'(4f6.2,a,2f6.2)') bet1,phi1,bet2,phi2, &
     &         ' WARNING: <Z>=',Kernel%zz(iq1,iq2,jproj,k1,k2) 
               if(abs(Nucl%nucleon(0)-Kernel%nn(iq1,iq2,jproj,k1,k2)).gt.0.1) &
     &         write(*,'(4f6.2,a,2f6.2)') bet1,phi1,bet2,phi2, &
     &         ' WARNING: <N>=',Kernel%nn(iq1,iq2,jproj,k1,k2) 
! .................... newly added for J^2
               read (lou,402) Kernel%J2(iq1,iq2,jproj,k1,k2)

               if((abs(bet1).gt.0.01 .and. abs(bet2).gt.0.01 ) .and. &
     &             abs(jproj*(jproj+1)-Kernel%J2(iq1,iq2,jproj,k1,k2)).gt.0.1) &
     &         write(*,'(4f6.2,a,i4,a,2f6.2)') bet1,phi1,bet2,phi2, &
     &         ' WARNING: should be',jproj*(jproj+1),&
     &                  ' but as ',Kernel%J2(iq1,iq2,jproj,k1,k2) 

              ! read (lou,402) Kernel%qp(iq1,iq2,jproj,k1,k2), Kernel%qcp(iq1,iq2,jproj,k1,k2)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,0),Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,1),Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
               read (lou,402) Kernel%qpred(iq1,iq2,jproj,k1,k2,2),Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
               read (lou,402) Kernel%q0p(iq1,iq2,jproj,k1,k2),Kernel%q0pc(iq1,iq2,jproj,k1,k2)

               !if(jproj.eq.2) print *,'Kernel%qpred <2||Q2||2>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.4) print *,'Kernel%qpred <4||Q2||4>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.6) print *,'Kernel%qpred <6||Q2||6>=',Kernel%qp(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.0) print *,'Kernel%qpred <2||Q2||0>=', Kernel%qpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.0) print *,'Kernel%qcpred<2||Q2||0>=',Kernel%qcpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.2) print *,'Kernel%qpred <4||Q2||2>=', Kernel%qpred(iq1,iq2,jproj,k1,k2)
               !if(jproj.eq.4) print *,'Kernel%qpred <6||Q2||4>=', Kernel%qpred(iq1,iq2,jproj,k1,k2)
!     ....................................... qcp, qcpred were set to zero for iq1=iq2 cases in PKC code
!          q0p    = <J(f),K1,q1 ||e*r^2|| J(i),K2,q2>
!          qp     = <J(f),K1,q1 ||e*Q_2|| J(i),K2,q2>
!          qcp    = <J,K1,q2   ||e*Q_2|| J,K2,q1>
!          qpred  = <J+2,K1,q1 ||e*Q_2|| J,K2,q2>
!          qcpred = <J+2,K1,q2 ||e*Q_2|| J,K2,q1>
!     ........................................ 
          if(iq1.ne.iq2) then
            !Kernel%qp(iq2,iq1,jproj,k1,k2)    = Kernel%qcp(iq1,iq2,jproj,k1,k2)
            Kernel%qpred(iq2,iq1,jproj,k1,k2,0) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,0)
            Kernel%qpred(iq2,iq1,jproj,k1,k2,1) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,1)
            Kernel%qpred(iq2,iq1,jproj,k1,k2,2) = Kernel%qcpred(iq1,iq2,jproj,k1,k2,2)
            Kernel%q0p(iq2,iq1,jproj,k1,k2)   = Kernel%q0pc(iq1,iq2,jproj,k1,k2)
            Kernel%nn(iq2,iq1,jproj,k1,k2)    = Kernel%nn(iq1,iq2,jproj,k1,k2)
            Kernel%zz(iq2,iq1,jproj,k1,k2)    = Kernel%zz(iq1,iq2,jproj,k1,k2)
          endif
!         enddo ! k2
!         enddo ! k1                           
      icheck = icheck + 1
  403 continue
!       write(*,*) 'icheck=',icheck 
  99  enddo ! jproj

      if(icheck.eq.0) then
      write(12,*) elem, 'does not exist ...'
      do jproj = jprmi,jprma,jdf
         k1  =  0
         k2  =  0
            Kernel%njkk(iq1,iq2,proj,k1,k2)   = 0.d0
            Kernel%hjkk(iq1,iq2,jproj,k1,k2)  = 0.d0
            !Kernel%qp(iq2,iq1,jproj,k1,k2)    = 0.d0
            !Kernel%qp(iq1,iq2,jproj,k1,k2)    = 0.d0
            Kernel%qpred(iq2,iq1,jproj,k1,k2,0:2) = 0.d0
            Kernel%qpred(iq1,iq2,jproj,k1,k2,0:2) = 0.d0
            Kernel%q0p(iq2,iq1,jproj,k1,k2)   = 0.d0
            Kernel%q0p(iq1,iq2,jproj,k1,k2)   = 0.d0
            Kernel%nn(iq1,iq2,jproj,k1,k2)= 0.d0
            Kernel%nn(iq2,iq1,jproj,k1,k2)= 0.d0
            Kernel%zz(iq2,iq1,jproj,k1,k2)= 0.d0
            Kernel%zz(iq1,iq2,jproj,k1,k2)= 0.d0
      enddo ! jproj
      endif
      return
      end
!.......................................................................
      subroutine Filename4Kernels(icase,iq1,iq2)
!.......................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-z)
!      implicit none
      INTEGER nuc1,nuc2,iq1,iq2,hw1,hw2
      character*1 sign1,sign2
      character*6 name0
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
!---------------------------  
              jphi    = PNP%NFOM 
              name81  = mod(jphi/10,10) + 48
              name82  = mod(jphi,10) + 48
              name91  = mod(AMP%NLEG_BET/10,10) + 48
              name92  = mod(AMP%NLEG_BET,10) + 48
              name_emax1 =mod(HO%emax/10,10) + 48
              name_emax2 =mod(HO%emax,10) + 48
              name_NOQ1  =mod(GCM%NOQ/10,10)+48
              name_NOQ2  =mod(GCM%NOQ,10)+48

              ncf11 = cf1+48
              ncf12 = mod(cf1*10,10.d0)  +48
              ncf13 = mod(cf1*100,10.d0) +48

              ncf21 = cf2+48
              ncf22 = mod(cf2*10,10.d0)  +48
              ncf23 = mod(cf2*100,10.d0) +48

              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48
              nuc1 = mod(Nucl%nucleon(2)/10,10) +48
              nuc2 = mod(Nucl%nucleon(2),10)    +48
!              if(mphi(1).lt.10)
!     &         name6  = '0'//char(jphi)
!---------------------------  
      if(Input%IsHFB .eq.2) then
       File%elem =name0//Nucl%nucnam//'/kern.VAP.'//char(name)//'D'&
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.elem'

      else 
      File%elem =name0//Nucl%nucnam//'/kern.'//char(name)//'D'         &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.elem'
      endif

      File%TD1B=name0//Nucl%nucnam//'/TD1B.'//Nucl%nucnam//char(nuc1)//char(nuc2) &
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
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

      File%Rho1B=name0//Nucl%nucnam//'/Rho1B.'//Nucl%nucnam//char(nuc1)//char(nuc2)&
     &           //'.'//char(name)//'D'    &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
!     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)          &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'
      File%Rho2B=name0//Nucl%nucnam//'/Rho2B.'//Nucl%nucnam//char(nuc1)//char(nuc2)&
     &           //'.'//char(name)//'D'     &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
!     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'
      File%Rho3B=name0//Nucl%nucnam//'/Rho3B.'//Nucl%nucnam//char(nuc1)//char(nuc2)&
     &           //'.'//char(name)//'D'     &
     &          //'.'//char(name81)//char(name82)            &
     &          //'.'//char(name91)//char(name92)//'.'       &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cValID)//'_'   &
!     &          //'_'//trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_' &
     &           //Input%cFlow//'_'                      &
     &          // trim(Input%vs4me3b)//'_'  &
     &          //sign1//char(name11)//char(name21)                 &
     &          //char(name31)//char(name41)//char(name51)//'.' &
     &          //char(name61)//char(name71)//char(namer1)    &
     &          //'cf'//char(ncf11)//char(ncf12)//char(ncf13)// &
     &          '_'                                             &
     &          //sign2//char(name12)//char(name22)//char(name32)&
     &          //char(name42)//char(name52)//'.'               &
     &          //char(name62)//char(name72)//char(namer2)      &
     &          //'cf'//char(ncf21)//char(ncf22)//char(ncf23)//'.dens'


      File%FF=name0//Nucl%nucnam//'/F4GS.'         &
!     &          //'eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //trim(Input%cIntID)//'_'   &
     &          //trim(Input%cValID)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &          //'NOQ'//char(name_NOQ1)//char(name_NOQ2)//'.dat'   

!--------------------- 
      return
      end
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
         if(.NOT. ALLOCATED(Const%hw_mesh))      ALLOCATE(Const%hw_mesh(1:NQ))
         do iq=1,NQ

            read(15,*) Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),&
                       Const%P00_mesh(iq), Const%hw_mesh(iq)

            !call Filename4wfs(Const%beta2t_mesh(iq),Const%gamma2t_mesh(iq),Const%P00_mesh(iq))
            !File%cwf(iq) = File%wf    ! name of files for wfs
!            print *,File%cwf(iq)
         enddo !
          Input%NGCM=NQ
!         print *,'Input%NGCM=',Input%NGCM
         if(Input%nq0f.gt.NQ .or. Input%nq1f.gt.NQ) &
     &   stop ' Error: Input%nq0f or Input%nq1f are larger than the total &
         &number of configurations!'
      end

        subroutine print_matrix(LDA,A_matrix)
        implicit none
        integer LDA,ii,jj
        real*8  A_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,f12.8)') jj,ii,A_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_zmatrix(LDA,ZA_matrix)
        implicit none
        integer     LDA,ii,jj
        complex*16  ZA_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,2e12.3)') jj,ii,ZA_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_array(LDA,V,istop)
        implicit none
        integer     LDA,ii,istop
        real*8 V(LDA)

        do ii=1,LDA
           write(911,'(i6,e12.3)') ii,V(ii)
        enddo
        if(istop .ne.0) stop

        END SUBROUTINE
!________________________________________________________________________
      subroutine Filename4wfs(betat,gammat,p00)
!.......................................................................
      USE VAPHFB_PAR
      implicit none
      character*1 sign1,sign2
      character*6 name0
      real*8      betat,gammat,p00,ab2c1
      integer     name_emax1,name_emax2,nuc1,nuc2,name1,name2,name3,&
     &            name4,name5,namep1,namep2,namep3,hw1,hw2
!---------------------------
              name0 = '../../'
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
              namep1 = p00+48
              namep2 = mod(p00*10,10.d0)  +48
              namep3 = mod(p00*100,10.d0) +48
              name_emax1 =mod(HO%emax/10,10) + 48
              name_emax2 =mod(HO%emax,10) + 48

              hw1= mod(Input%ihwHO/10,10)+48
              hw2= mod(Input%ihwHO,10)+48
!---------------------------  
       if(Input%IsHFB.eq.0) then

      File%wf =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_HF_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)//'.wf' ! 


      else if(Input%IsHFB.eq.2) then
      File%wf =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_VAP_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)//'.wf' !      

      else if(Input%IsHFB.eq.1) then
      File%wf =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)    &
     &           //'_HFB_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &          //'np'//char(namep1)//char(namep2)//'.wf' !      
      endif

!---------------------------  
      if(Input%IsHFB.eq.2) then
      File%out =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_VAP_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &          //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      

      else if(Input%IsHFB.eq.1) then
      File%out =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_HFB_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &        //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      

      else if(Input%IsHFB.eq.0) then
      File%out =name0//Nucl%nucnam//'/'//Nucl%nucnam//char(nuc1)//char(nuc2)   &
     &           //'_HF_'//Input%cIntID       &
     &          //'_eMax'//char(name_emax1)//char(name_emax2)//'_'   &
     &          //'hwHO'//char(hw1)//char(hw2)//'_'   &
     &           //Input%cFlow//'_'                      &
     &           //'beta'//sign1//char(name1)                       &
     &          //char(name2)//char(name3)                           &
     &        //'gam'//char(name4)//char(name5)                    &
     &       //'np'//char(namep1)//char(namep2)//'.out' !      
      endif
     
      return
      end





!=======================================================================
      subroutine nucleus(is,npro,te)

!=======================================================================
!
!     is = 1 determines the symbol for a given proton number npro
!          2 determines the proton number for a given symbol te
!
!-----------------------------------------------------------------------
!
      PARAMETER (MAXZ=140)
!
      CHARACTER TE*2,T*(2*MAXZ+2)
!
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
!
! ... Rf is called also as Ku (kurchatovium)
! ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
! ... Ha is called also as Db (Dubnium)
!
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong'
         te = t(2*npro+1:2*npro+2)
         return
      else
!
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
               if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
!
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
!
      stop
      END
!
