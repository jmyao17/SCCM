       subroutine Calc_TD1B(iq1,iq2,lpr)
       USE VAPHFB_Par
       implicit none
       logical lpr
       integer n12,iq1,iq2
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,nn,np
       real*8 TD,TD1B
       real*8 TDc,TD1Bc
       dimension TD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)
       dimension TD1Bc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)

       integer ji,jf,li,lf,Lam,tnlj1,tnlj2

!      .... read density
       call read_TD1B(33,TD1B,TD1Bc)

!      .... sum density over q1,q2

       ! ki=kf=1 (first state for a given J)
       do Ji=0,GCM%Jmax
       do li=1,GCM%kmax
       do Jf=0,GCM%Jmax
       do lf=1,GCM%kmax
          !n12 = 1
          n12 = 1 + iq1/iq2  ! 1 or 2
       do Lam=0,2,2
       do tnlj1=1,FSPB%tnlj_fmax
       do tnlj2=1,FSPB%tnlj_fmax
          GCMDens%TD1BSum(tnlj1,tnlj2,Ji,Jf,li,lf,Lam) = &
      &   GCMDens%TD1BSum(tnlj1,tnlj2,Ji,Jf,li,lf,Lam) &
      &                           + ( GCM%FqJk(iq1,Jf,lf)*GCM%FqJk(iq2,Ji,li)& 
      &                                *TD1B(tnlj1,tnlj2,Ji,Jf,Lam)  &
      &                           + GCM%FqJk(iq2,Jf,lf)*GCM%FqJk(iq1,Ji,li) & 
      &                              *TD1Bc(tnlj1,tnlj2,Ji,Jf,Lam))/n12

       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
!     ............
       end

       subroutine Write_TD1B(lpr)
       USE VAPHFB_Par
       implicit none
       logical lpr
       integer e0,it,icc,n1,n2,L1,L2,J1,J2,LJ,nn,np
       real*8 TD,TD1B
       real*8 TDc,TD1Bc
       integer ji,jf,li,lf,Lam,tnlj1,tnlj2
       integer jt1,lj1,jt2,lj2,jj1,jj2
        integer    SPB_twoj_lj,SPB_l_lj

       open(912,file='GCM_TD1BSum.dat',status='unknown')
       write(912,*) ' t   J_i   l_i   J_f   l_f',&
                   '   Lam  n1   lj1  n2  lj2', &
                   '   TD   TD/sqrt(2j_1+1)'
       if(.not. lpr) return
       do Ji=0,GCM%Jmax
       do li=1,GCM%kmax
       do Jf=0,GCM%Jmax
       do lf=1,GCM%kmax
       do Lam=0,2,2
       do tnlj1=1,FSPB%tnlj_fmax

           jt1 = FSPB%t(tnlj1) ! n(0), p(1)
           n1 = FSPB%n(tnlj1) ! starting from 1 
           lj1 = FSPB%lj(tnlj1)
           jj1 = SPB_twoj_lj(lj1)    ! doubled
       do tnlj2=1,FSPB%tnlj_fmax
           jt2 = FSPB%t(tnlj2) ! n(0), p(1)
           n2 = FSPB%n(tnlj2) ! starting from 1 
           lj2 = FSPB%lj(tnlj2)
           jj2 = SPB_twoj_lj(lj2)    ! doubled
           TD  = GCMDens%TD1BSum(tnlj1,tnlj2,Ji,Jf,li,lf,Lam)
           if(abs(TD).gt.1.d-5) &
           write(912,'(i3,4i6,i3,4i6,2f15.8)') &
           jt1,Ji,li,Jf,lf,Lam,n1,lj1,n2,lj2,TD,TD/sqrt(jj1+1.0)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       print *, ' -> transition density is written into GCM_TD1BSum.dat...'
!     ............
       end

       subroutine read_TD1B(ie,TD1B,TD1Bc)
       USE VAPHFB_PAR
       implicit none
!      ...............................
!      unnormalized transition density
!      ...............................
       integer ie
       real*8 TD,TD1B
       real*8 TDc,TD1Bc
       dimension TD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  !
       dimension TD1Bc(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  !
       integer Ji,Jf,Lam,mu,tnlj1,tnlj2

       TD1B(:,:,:,:,:)  = 0.d0 
       TD1Bc(:,:,:,:,:) = 0.d0 
       open(ie,file=File%TD1B,status='old')
       ! skip the first five lines
       read(ie,*) 
       read(ie,*) 
       read(ie,*) 
       read(ie,*) 
       read(ie,*) 
  7    read(ie,10,end=8) Ji,Jf,Lam,tnlj1,tnlj2,TD,TDc
       TD1B(tnlj1,tnlj2,Ji,Jf,Lam)  = TD 
       TD1Bc(tnlj1,tnlj2,Ji,Jf,Lam) = TDc 
       goto 7
  8   continue
      close(ie) 
      return
10     format(5i8,2f15.8)
       end
