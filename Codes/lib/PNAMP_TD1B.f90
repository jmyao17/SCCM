!      ..................................................................      
!       computing the reduced transition density: Rho_AB (Lam: Ji -> Jf)
!      ..................................................................      
        subroutine ZTD1B_Angles(iangle,ZRO,Znorm,ZTD1B,ZTD1Bc,NLEV)
        USE VAPHFB_PAR
        implicit none
        integer   mitf,iangle,NLEV 
        integer   k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 ZRO,ZTD1B,ZTD1Bc,Znorm
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION ZTD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)  
        DIMENSION ZTD1Bc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)  
        integer Lam,mu,tnlj1,tnlj2
        real*8 cg1

!     ............. print out the one-body density
        do k1=1,NLEV
        do k2=1,NLEV
           lj1 = tnljm%lj(k1)
           lj2 = tnljm%lj(k2)
           twoj1  = tnljm%twoj(k1)
           twoj2  = tnljm%twoj(k2)
           twom1  = tnljm%twom(k1)
           twom2  = tnljm%twom(k2)
           it1  = mitf(tnljm%t(k1))        ! tnljm%t(k1)=-1 for p
           it2  = mitf(tnljm%t(k2))
           n1   = tnljm%n(k1)              ! n=1,2,3,...
           n2   = tnljm%n(k2)
           
           tnlj1=FSPB%tnlj(lj1,n1,it1)     ! it1=0(n),1(p)
           tnlj2=FSPB%tnlj(lj2,n2,it2)     ! it1=0(n),1(p)

            ! .............................................................
            ! znorm: <F| R(Omega, varphi_n, varphi_p)| I>
            ! ZRO(k2,k1): <F| a^+_k1  a_k2| I (Omega, varphi_n, varphi_p)> 
            ! ............................................................

            do Lam=0,2,2
               do mu=-Lam,Lam
                  ! j1+0.5, m1+0.5, j2+0.5, -m2+0.5, J12, M12 
                  cg1 = CG_Save((twoj1+1)/2,(twom1+1)/2,(twoj2+1)/2,(-twom2+1)/2,Lam,mu)

                  ! < q_f | ([a^+_A \tilde a_B]_Lam,mu ) R(Omega)|q_i>
                  ZTD1B(tnlj1,tnlj2,Lam,mu)= ZTD1B(tnlj1,tnlj2,Lam,mu) &
     &            + ZRO(k2,k1) *cg1*znorm*iv(abs(twoj2-twom2)/2)/PNP%NFOM**2  

                  ! < q_i | ([a^+_A \tilde a_B]_Lam,mu )^+ R(Omega)|q_f>*
                  ZTD1Bc(tnlj1,tnlj2,Lam,mu)= ZTD1Bc(tnlj1,tnlj2,Lam,mu) &
     &            + ZRO(k1,k2)*cg1*znorm*iv(abs(twoj2-twom2)/2)/PNP%NFOM**2 

               enddo
            enddo
!      ......................
        enddo
        enddo
        return

        ! ..... for test
       do Lam=0,2,2
       do mu=-Lam,Lam
       do tnlj1=1,FSPB%tnlj_fmax
       do tnlj2=1,FSPB%tnlj_fmax
          if(abs(ZTD1B(tnlj1,tnlj2,Lam,mu)).gt.1.d-4) &
          write(*,'(2i6,2i4,2f12.8)') tnlj1,tnlj2,Lam,mu,dreal(ZTD1B(tnlj1,tnlj2,Lam,mu)), &
               dreal(ZTD1Bc(tnlj1,tnlj2,Lam,mu))
       enddo
       enddo
       enddo
       enddo
        stop
       end


        subroutine ZTD1B_Integral(alp,bet,gam,ZTD,ZTDc,zff,ZTDF,ZTDFc)
        USE VAPHFB_PAR
        implicit none
        integer   iangle
        REAL(DP) CG
        real*8    dwignerI_gen
        integer   J,k1,k2,it1,it2,n1,n2,twom1,twom2,twoj1,twoj2,lj1,lj2
        complex*16 zap,zgp,zbp,zan,zgn,zbn,zff,ZTDF,ZTD,ZTDc,ZTDFc
        real*8  cg1,bet,alp,gam
        dimension ZTD(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)  
        dimension ZTDc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:2,-2:2)  
        dimension ZTDF(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  ! Ki=Kf=0
        dimension ZTDFc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  ! Ki=Kf=0
        integer Ji,Jf,Lam,nu,mu,tnlj1,tnlj2
        integer Ki,Kf


      
        Ki=0
        Kf=0
       do Lam=0,2,2
       do mu=-Lam,Lam

          zap=cdexp(Zimag*(Kf-mu)*alp)
          zgp=cdexp(Zimag*Ki*gam)
       do Ji=0,JJmax
          zbp = zone*dwignerI_gen('i',Ji,Kf-mu,Ki,bet)
       do Jf=0,JJmax
          cg1 = CG(Ji*2,2*(Kf-mu),Lam*2,mu*2,Jf*2,Kf*2) 

       do tnlj1=1,FSPB%tnlj_fmax
       do tnlj2=1,FSPB%tnlj_fmax

          ! zff: measure for integration
         ZTDF(tnlj1,tnlj2,Ji,Jf,Lam)= ZTDF(tnlj1,tnlj2,Ji,Jf,Lam) &
     &                        + iv(2*Lam)*sqrt(2*Jf+1.0)*cg1              &
     &                         *zap*zbp*zgp*weightJ(Ji)*zff    &
                               *ZTD(tnlj1,tnlj2,Lam,mu) 

       enddo
       enddo
       enddo
       enddo
       enddo
       enddo


!       ....... exchage qi and qf

       do Lam=0,2,2
       do mu=-Lam,Lam ! nu
       do nu=-Lam,Lam ! nu' 

          zan=cdexp(-Zimag*(Ki+nu)*alp)
          zgn=cdexp(-Zimag*(Kf-2*mu)*gam)
       do Ji=0,JJmax
          zbn = zone*dwignerI_gen('i',Ji,Ki,Kf-mu,bet)*dwignerI_gen('i',Lam,nu,-mu,bet)
       do Jf=0,JJmax
          cg1 = CG(Ji*2,2*(Kf-mu),Lam*2,mu*2,Jf*2,Kf*2)

       do tnlj1=1,FSPB%tnlj_fmax
       do tnlj2=1,FSPB%tnlj_fmax

          ! zff: measure for integration
         ZTDFc(tnlj1,tnlj2,Ji,Jf,Lam)= ZTDFc(tnlj1,tnlj2,Ji,Jf,Lam) &
     &                        + iv(Lam-mu)*sqrt(2*Jf+1.0)*cg1   &
     &                         *zan*zbn*zgn*weightJ(Ji)*zff    &
                               *dconjg(ZTDc(tnlj1,tnlj2,Lam,nu))

       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

       end

       subroutine Write_TD1B(ie,ZTD1B,ZTD1Bc,znorm) 
       USE VAPHFB_PAR
       implicit none
!      ...............................
!      unnormalized transition density
!      ...............................
       integer ie      
       complex*16 zTD,zTDc,znorm
       complex*16 ZTD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  ! Ki=Kf=0
       complex*16 ZTD1Bc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  ! Ki=Kf=0
       integer Ji,Jf,Lam,mu,tnlj1,tnlj2
       write(ie,*)  '---------------------------------------------'
       write(ie,*) '      TD(qf,qi): <Jf Kf qf|| T_L || Ji Ki qi> '
       write(ie,*) '      TD(qi,qf): <Jf Kf qi|| T_L || Ji Ki qf>'
       write(ie,*)  '---------------------------------------------'
       write(ie,*) '      Ji     Jf     Lam     tnlj1    tnlj2      TD(qf,qi)   TD(qi,qf)'
       do Ji=0,JJmax
       do Jf=0,JJmax
       do Lam=0,2,2
       do tnlj1=1,FSPB%tnlj_fmax
       do tnlj2=1,FSPB%tnlj_fmax
          zTD  = ZTD1B(tnlj1,tnlj2,Ji,Jf,Lam)  ! /znorm
          zTDc = ZTD1Bc(tnlj1,tnlj2,Ji,Jf,Lam) ! /znorm
         if(abs(zTD).gt.1.d-9 .or. abs(zTDc).gt.1.d-9) write(ie,10) Ji,Jf,Lam,tnlj1,tnlj2,dreal(zTD),dreal(zTDc)

       enddo
       enddo
       enddo
       enddo
       enddo
10     format(5i8,2f15.8)
       end

