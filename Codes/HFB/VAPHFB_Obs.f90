!    .......................................
        subroutine Calc_Obs_From_UV(rank,ie,epsi,U,V,RO_0,akapa10_0,akapa01_0, &
     &             akin,AZ_P,AN_P,Q2_2t,Q2_1t,Q20t,Q21t,Q22t,     &
     &             AJX_ME,AJX2_ME,AJY_ME,AJY2_ME,AJZ_ME,AJZ2_ME,  &
     &             P_1m1_00_me,P_10_00_me,P_1p1_00_me,  &
     &             P_00_1m1_me,P_00_10_me,P_00_1p1_me,NLEV,NLEV2)
!    .......................................
        USE VAPHFB_PAR


        implicit real*8 (a-h,o-z)

        integer :: rank
        !DENSITIES AND FIELDS
        DIMENSION RO_0(NLEV,NLEV)       ! density matrix
        DIMENSION akapa10_0(NLEV,NLEV)        ! pairing tensor 10 (cc)
        DIMENSION akapa01_0(NLEV,NLEV)        ! pairing tensor 01 (c+c+)
        !WAVE FUNCTIONS
        DIMENSION U(NLEV,NLEV)
        DIMENSION V(NLEV,NLEV)

        !OTHER QUANTITIES
        DIMENSION akin(NLEV,NLEV)
!        DIMENSION AZme(NLEV2),ANme(NLEV2)
        DIMENSION Q20t(NLEV2)
        DIMENSION Q22t(NLEV2),Q21t(NLEV2)
        DIMENSION Q2_2t(NLEV2),Q2_1t(NLEV2)
        DIMENSION AJX_ME(NLEV2),AJX2_ME(NLEV2)
        DIMENSION AJY_ME(NLEV2),AJY2_ME(NLEV2)
        DIMENSION AJZ_ME(NLEV2),AJZ2_ME(NLEV2)
        DIMENSION P_00_1m1_me(NLEV2),P_00_1p1_me(NLEV2)
        DIMENSION P_1m1_00_me(NLEV2),P_1p1_00_me(NLEV2)
        DIMENSION P_00_10_me(NLEV2),P_10_00_me(NLEV2)
        Real*8    Etot,epsi

1     format('**************************************')
2     format('Precise:',f12.8)
      write(ie,1) 
      write(ie,'(a)') '  Summary of Results ' 
      write(ie,1) 
      write(ie,2) epsi
      write(ie,1) 

      write(*,*) ' computing energy with PAV ...'
!       Energy 


        call NZPROJ_PAV(U,V,H%ME1BM,AZ_P,AN_P,            &
     &  EKin_P,Ekin_N,                                 &
     &  EHF_PP,EHF_PN,EHF_NP,EHF_NN,                   &
     &  EPa_PP,EPa_PN,EPa_NP,EPa_NN,                   &
     &  EHFB_P,EHFB_N,PNP%NFOM,NLEV)

        !energy !CHANGE FOR PNP-energies.

        write(ie,'(A20,3A15)') ' ','proton','neutron',      &
     & 'total'
        write(ie,'(A20,3A15)') ' ','------','------','------'
        write(ie,*) ' '

        write(ie,'(A20,3F15.6)') 'Kinetic Energy',Ekin_P,Ekin_N,  &
     &  Ekin_P+Ekin_N
        write(ie,'(A20,3F15.6)') 'PNP Energy',EHFB_P,EHFB_N,EHFB_P+EHFB_N


        write(ie,*) ' '
        write(ie,'(A20,5A15)') ' ','p-p','n-n','p-n','n-p',  &
     &  'total'
        write(ie,'(A20,5A15)') ' ','------','------','------',  &
     &  '------','------'
        write(ie,*) ' '

        write(ie,'(A20,5F15.6)') 'Kin Energy',Ekin_P,Ekin_N,0.d0,0.d0, &
     &  Ekin_P+Ekin_N

        write(ie,'(A20,5F15.6)') 'HF Energy',   &
     &  EHF_PP,EHF_NN,EHF_PN,EHF_NP,           &
     &  EHF_PP+EHF_NN+EHF_PN+EHF_NP

        write(ie,'(A20,5F15.6)') 'Pair Energy',  &
     &  EPa_PP,EPa_NN,EPa_PN,EPa_NP,            &
     &  EPa_PP+EPa_NN+EPa_PN+EPa_NP


!     ......... pairing for neutrons
          if(abs(EPa_NN).gt.0.0001) then 
            HFB%ide(0) =1
          else
            HFB%ide(0) =0
          endif

!     ......... pairing for protons 
          if(abs(EPa_PP).gt.0.0001) then 
            HFB%ide(1) =1
          else
            HFB%ide(1) =0
          endif

!     ......... pairing for n-p
          if(abs(EPa_PN).gt.0.00001) then
            print *,' np isoscalar pairing exists ....'
            HFB%ide(0) =1
            HFB%ide(1) =1
          endif


        write(ie,'(A20,a15,a15,a15,a15,f15.6)') &
     &  'Zero-Point-Energy','', '', '', '', H%E0

       write(ie,'(A100)') '---------------------'
       write(ie,'(A20,a15,a15,a15,a15,f15.6)') 'Sum','', '', '', '', Ekin_P+Ekin_N+H%E0 &
     &     +EPa_PP+EPa_NN+EPa_PN+EPa_NP+EHF_PP+EHF_NN+EHF_PN+EHF_NP

        Etot = Ekin_P+Ekin_N+H%E0                                    &
     &     +EPa_PP+EPa_NN+EPa_PN+EPa_NP+EHF_PP+EHF_NN+EHF_PN+EHF_NP
!       ............
!       deformation
!       ............
        write(ie,*) ' '
        write(ie,*) 'QUADRUPOLE DEFORMATIONS'
        write(ie,*) ' '

!      ........ deformation parameter
        call Q_0(Q20t,ro_0,Q20p,Q20n,NLEV)
        write(ie,'(A20,3F15.6)') 'Q20',Q20p,Q20n,Q20p+Q20n
        Q20nbeta = Q20n
        Q20pbeta = Q20p
        Q20tbeta = Q20p+Q20n
        call Q_0(Q2_2t,ro_0,Q2_2p,Q2_2n,NLEV)
        write(ie,'(A20,3F15.6)') 'Q2_2',Q2_2p,Q2_2n,Q2_2p+Q2_2n
        Q2_2nbeta = Q2_2n
        Q2_2pbeta = Q2_2p
        Q2_2tbeta = Q2_2p+Q2_2n
        call Q_0(Q22t,ro_0,Q22p,Q22n,NLEV)
        write(ie,'(A20,3F15.6)') 'Q22',Q22p,Q22n,Q22p+Q22n
        Q22nbeta = Q22n
        Q22pbeta = Q22p
        Q22tbeta = Q22p+Q22n
!     ............. determination of beta2 from Q20 and Q22
        bsignn = Q20nbeta/abs(Q20nbeta)
        bsignp = Q20pbeta/abs(Q20pbeta)
        bsignt = Q20tbeta/abs(Q20tbeta)
        beta20n=bsignn*Const%Q2BN * dsqrt(Q20nbeta**2+Q22nbeta**2+Q2_2nbeta**2)
        beta20p=bsignp*Const%Q2BP * dsqrt(Q20pbeta**2+Q22pbeta**2+Q2_2pbeta**2)
        beta20t=bsignt*Const%Q2BA * dsqrt(Q20tbeta**2+Q22tbeta**2+Q2_2tbeta**2)

        write(ie,'(A20,3F15.6)') 'beta20',beta20p,beta20n,beta20t
        gam2p  = 180.d0*atan(sqrt(2.d0)*Q22p/Q20pbeta)/pi
        gam2n  = 180.d0*atan(sqrt(2.d0)*Q22n/Q20nbeta)/pi
        gam2t  = 180.d0*atan(sqrt(2.d0)*Q22tbeta/Q20tbeta)/pi
        write(ie,'(A20,3F15.6)') 'gamma22',gam2p,gam2n,gam2t

        call Q_0(Q2_1t,ro_0,Q_0p,Q_0n,NLEV)
         write(ie,'(A20,3F15.6)') 'Q2_1',Q_0p,Q_0n,Q_0p+Q_0n
        call Q_0(Q21t,ro_0,Q_0p,Q_0n,NLEV)
         write(ie,'(A20,3F15.6)') 'Q21',Q_0p,Q_0n,Q_0p+Q_0n

        write(ie,*) 'rms radii'
        call Q_0(cME1B%r2,ro_0,r2p,r2n,NLEV)
        write(ie,'(A20,3F15.6)') 'rms radii',dsqrt(r2p/Input%nprot),dsqrt(r2n/Input%nneut),&
                  &          dsqrt((r2p+r2n)/(Input%nprot+Input%nneut))

        !Jx,Jy,Jz
        write(ie,*) ' '
        write(ie,*) 'AXIALITY'
        write(ie,*) ' '

!     ................ Jx
        call Q_0(AJX_ME,ro_0,Q_0p,Q_0n,NLEV)
        call F2_0(RO_0,AKAPA10_0,AKAPA01_0,AJX_ME,AJX2_ME,AJX2_MV,NLEV)
         write(ie,'(A20,3F15.6)') 'JX JX2 DJX',  &
     &   Q_0p+Q_0n,AJX2_MV,AJX2_MV-(Q_0p+Q_0n)**2

        call Q_0(AJY_ME,ro_0,Q_0p,Q_0n,NLEV) ! check imaginary Jy!
        call F2_0(RO_0,AKAPA10_0,AKAPA01_0,      &
     &  AJY_ME,AJY2_ME,AJY2_MV,NLEV)
        
         write(ie,'(A20,3F15.6)') 'JY JY2 DJY',         &
     &   -(Q_0p+Q_0n),-AJY2_MV,-AJY2_MV-(Q_0p+Q_0n)**2


        call Q_0(AJZ_ME,ro_0,Q_0p,Q_0n,NLEV)
        call F2_0(RO_0,AKAPA10_0,AKAPA01_0,AJZ_ME,AJZ2_ME,AJZ2_MV,NLEV)
         write(ie,'(A20,3F15.6)') 'JZ JZ2 DJZ',          &
     &    Q_0p+Q_0n,AJZ2_MV,AJZ2_MV-(Q_0p+Q_0n)**2


!      ........ deformation parameter
        write(ie,*) 'Octuple Deformation'
        call Q_0(cME1B%Q30t,ro_0,Q30p,Q30n,NLEV)
        write(ie,'(A20,3F15.6)') 'Q30',Q30p,Q30n,Q30p+Q30n


        write(ie,*) ' '
        write(ie,*) '///////////////////////'
        write(ie,*) '         PAIRS'
        write(ie,*) '/////////////////////// '
        write(ie,*) ' '
        call Qpair_0(P_00_1m1_me,Akapa10_0,Akapa01_0,Q_T00_J1m1,NLEV)
        call Qpair_0(P_00_10_me,Akapa10_0,Akapa01_0,Q_T00_J10,NLEV)
        call Qpair_0(P_00_1p1_me,Akapa10_0,Akapa01_0,Q_T00_J1p1,NLEV)

        write(ie,'(A20,3A15)') ' ','MJ=-1','MJ=0','MJ=+1'
        write(ie,'(A20,3A15)') ' ','------','------','------'
        write(ie,*) ' '
         write(ie,'(A20,3F15.6)') 'T=0 J=1',                   &
     &   dabs(Q_T00_J1m1),dabs(Q_T00_J10),dabs(Q_T00_J1p1)

        !T=1 ; J=0
        call Qpair_0(P_1m1_00_me,Akapa10_0,Akapa01_0,Q_T1m1_J00,NLEV)
        call Qpair_0(P_10_00_me,Akapa10_0,Akapa01_0,Q_T10_J00,NLEV)
        call Qpair_0(P_1p1_00_me,Akapa10_0,Akapa01_0,Q_T1p1_J00,NLEV)
        write(ie,*) ' '
        write(ie,'(A20,3A15)') ' ','MT=-1 PP','MT=0 PN','MT=+1 NN'
        write(ie,'(A20,3A15)') ' ','------','------','------'
        write(ie,*) ' '
        write(ie,'(A20,3F15.6)') 'T=1 J=0',                &
     &  dabs(Q_T1m1_J00),dabs(Q_T10_J00),dabs(Q_T1p1_J00)




        write(900+rank,'(f10.6,4f7.3,3f9.3,3f9.3,2f10.3)') &
     &  epsi,shared%cf,beta20n,beta20p,beta20t, &
     &  gam2n,gam2p,gam2t,dabs(Q_T00_J10),EPa_NN,EPa_PP,EPa_PN+EPa_NP,Etot


        end subroutine


        SUBROUTINE canon_basis(io,ro,gamma,delta,ham,NDIM)
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        DIMENSION gamma(NDIM,NDIM),delta(NDIM,NDIM),ham(NDIM,NDIM)
        DIMENSION ro(NDIM,NDIM)
        DIMENSION AUX9(NDIM,NDIM),AUX10(NDIM,NDIM)
        DIMENSION AUX11(NDIM,NDIM),AUX12(NDIM,NDIM)
        DIMENSION occup(NDIM,NDIM)
        DIMENSION Ocu_rho(NDIM)
        DIMENSION spe_ham(NDIM)
        DIMENSION WORKrho(3*NDIM-1)
        DIMENSION WORKham(3*NDIM-1)
        DIMENSION roaux(NDIM,NDIM)
        DIMENSION hamaux(NDIM,NDIM)
        real*8 hh(NDIM,NDIM),ee(NDIM),ez(NDIM)
        real*8 dkkbar(NDIM,NDIM),de(0:2)
!       Fermi energy

        write(io,'(A20,2F15.6)') 'Fermi Ener(Lambda):',HFB%EFermi(0),HFB%EFermi(1) 
        write(io,*) '    '
        de(0:2)=zero
        do ii=1,NDIM
         do jj=1,NDIM
          hamaux(jj,ii) = zero 
          roaux(jj,ii)  =-ro(jj,ii)        ! minus sign is purely for ordering the levels
          hamaux(jj,ii) = H%ME1BM(jj,ii)+gamma(jj,ii)
          hh(jj,ii)     = H%ME1BM(jj,ii)+gamma(jj,ii)
          de(2) = de(2) + delta(ii,jj)*ro(jj,ii)
         end do
        end do



!        call sdiag(NDIM,NDIM,roaux,Ocu_rho,roaux,ez,1)

!  DSYEV: computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A        
!  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )            
!  ----------------------------
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!  A       (input/output) DOUBLE PRECISION array
!  W       (output) eigenvalues
!       ........ diagonal density matrix rho -> U_rho
        call dsyev('V','U',NDIM,roaux,NDIM,Ocu_rho,WORKrho,3*NDIM-1,INFO2)

!       ........ diagonal h matrix          -> U_h
        call dsyev('V','U',NDIM,hamaux,NDIM,spe_ham,WORKham,3*NDIM-1,INFO3)

! SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!> DGEMM  performs one of the matrix-matrix operations
!    C := alpha*op( A )*op( B ) + beta*C,
!   32 *> where  op( X ) is one of
!   33 *>
!   34 *>    op( X ) = X   or   op( X ) = X**T,
!   35 *>
!   36 *> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!   37 *> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!   40 *  Arguments:
!   41 *  ==========
!   45 *>          TRANSA is CHARACTER*1
!   46 *>           On entry, TRANSA specifies the form of op( A ) to be used in
!   47 *>           the matrix multiplication as follows:
!   48 *>
!   49 *>              TRANSA = 'N' or 'n',  op( A ) = A.
!   51 *>              TRANSA = 'T' or 't',  op( A ) = A**T.
!   53 *>              TRANSA = 'C' or 'c',  op( A ) = A**T.


!    Using  U_rho to transform h
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,   &
     &  roaux,NDIM,ham,NDIM,zero,AUX9,NDIM)

        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX9,NDIM,roaux,NDIM,zero,AUX10,NDIM)

!    Using  U_h to transform rho 
!        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
!     &  hamaux,NDIM,ro,NDIM,zero,AUX11,NDIM)
!     &  roaux,NDIM,ro,NDIM,zero,AUX11,NDIM)
!        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
!     &  AUX11,NDIM,roaux,NDIM,zero,occup,NDIM)
!     &  AUX11,NDIM,hamaux,NDIM,zero,occup,NDIM)


!        do ii=1,NDIM
!            write(io,'(i5,F15.7)') ii,-Ocu_rho(ii)
!        enddo

        write(io,*) '-----------------------------------------'
        write(io,*) ' occ. from diagonalization of V*V^T, and      '
        write(io,*) ' the same transf. is applied onto h.   '
        write(io,*) '-----------------------------------------'
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  roaux,NDIM,t3me,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,roaux,NDIM,zero,AUX12,NDIM)

!   .......................... transform D_m1m2 to D_k1k2
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  roaux,NDIM,delta,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,roaux,NDIM,zero,dkkbar,NDIM)
!   ............................... print out s.p. levels with quantum number
 103    format(i3,'  ',a20,'   ',f5.2,'  ',f8.3, '  ',3f12.7)
 104    format(//,' single-particle properties of neutrons ',/,1x,33(1h-))
         write(io,104) 
         write(io,'(a)') ' idx quantum numbers        coef       ee          vv  &
         &        tz       delta(k,kbar)'
        do kk=1,HO%NLEV
           cmax = zero     
           imax = 0
           if(abs(AUX12(kk,kk)-1).gt.0.5) cycle
           do mm=1,HO%NLEV
              if(abs(roaux(mm,kk)).gt. abs(cmax)) then
                 cmax = roaux(mm,kk)    
                 imax = mm    
              endif
           enddo !
           if(abs(Ocu_rho(kk)).gt.1.d-5) then
             write(io,103) kk,HO%tts(imax),cmax,AUX10(kk,kk),-Ocu_rho(kk),AUX12(kk,kk),dkkbar(kk,kk-iv(kk))
           endif
             de(0) = de(0)+abs(Ocu_rho(kk)*dkkbar(kk,kk-iv(kk)))
        enddo !


        de(0)=0.d0
        de(1)=0.d0
        de(2)=0.d0

 105    format(//,' single-particle properties of protons ',/,1x,33(1h-))
         write(io,105)
         write(io,'(a)') ' idx quantum numbers  coef     ee       vv'
        do kk=1,HO%NLEV
           cmax = zero
           imax = 0
           de(2) = de(2)+abs(Ocu_rho(kk)*dkkbar(kk,kk-iv(kk)))
!          ...... loop over only proton
           if(abs(AUX12(kk,kk)+1).gt.0.5) cycle
           do mm=1,HO%NLEV
              if(abs(roaux(mm,kk)).gt. abs(cmax)) then
                cmax = roaux(mm,kk)
                imax = mm
              endif
           enddo !
           if(abs(Ocu_rho(kk)).gt.1.d-5) then
             write(io,103) kk,HO%tts(imax),cmax,AUX10(kk,kk),-Ocu_rho(kk),AUX12(kk,kk),dkkbar(kk,kk-iv(kk))
           endif
             de(1) = de(1)+abs(Ocu_rho(kk)*dkkbar(kk,kk-iv(kk)))
        enddo !
        write(io,*) ' ' 
        write(io,*) '--------------------------------------------------------------'
        write(io,'(A30,3F10.6)') 'Average Pairing Gap(n/p/tot):',de(0)/Nucl%nucleon(0),&
       &    de(1)/Nucl%nucleon(1),de(2)/Nucl%nucleon(2) 
        write(io,*) '--------------------------------------------------------------'
       

! ........... method 2
 106    format(//,' single-particle properties of neutrons/protons ',/,1x,33(1h-))
        call sdiag(NDIM,NDIM,hh,ee,hh,ez,1)

!    Using  U_h to transform t3 
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  hh,NDIM,t3me,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,hh,NDIM,zero,AUX12,NDIM)

!    Using  U_h to transform rho 
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,    &
     &  hh,NDIM,ro,NDIM,zero,AUX11,NDIM)
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,    &
     &  AUX11,NDIM,hh,NDIM,zero,occup,NDIM)

        write(io,106)

        do it=-1,1,2
        if(it.eq.-1) write(io,*) 'Proton:'
        if(it.eq. 1) write(io,*) 'Neutron:'
        do kk=1,NDIM
           cmax = zero
           imax = 0
        do mm=1,NDIM
           if(abs(hh(mm,kk)).gt. abs(cmax)) then
             cmax = hh(mm,kk)
             imax = mm
           endif
        enddo ! mm
          if(ee(kk).lt.50.d0 .and. abs(AUX12(kk,kk)-it).lt.0.1) then
             write(io,103) kk,HO%tts(imax),cmax,ee(kk),AUX12(kk,kk),occup(kk,kk)
          endif 
        enddo ! kk
        enddo
!       call checkwf(HFB%U0,HFB%V0,HO%NLEV)
!        print *, ' Single-particle energy levels have been written into file'
        return


        END SUBROUTINE


