!______________________________________________________________________________
      subroutine trans_by_td(jmn,jmx,jdf)

!..............................................................................
!     calculate reduced transition rates using transition density             .
!                                                                             .
!     iqk = iq+k1m/2*maxmp
!     k1m = k1-iki0, where iki0=0 (even J) or 2 (odd J) 
!     iq = number of points in q space                                        .
!     ik = number of points in k space                                        .
!     k  = number of eigenvalues of the norm kernel                           . 
!                                                                             .
!.............................................................................. 
!     EJKq(k,J)    : energy for state J_k                                     .
!     FJKq(iqk,k,J): wave function for state J_k                              .
!..............................................................................
      USE VAPHFB_Par
      implicit real*8(a-h,o-z)
      
!..............................................................................
   20 format(2i2,3f12.6) 
  203 format (4i5,1f10.2,1f12.5,1f8.3,2f10.3)
  202 format (' i1 Ji i2 Jf   Ei    Ef    DeltaE',               &
     &        '  E2:(J_i->J_f)(J_f -> J_i) |Q(t)|',              &
     &        ' |beta(t)| <Ji||E2||Jf>   <er^2>    rho^2(E0)',/,  &
     &        '               MeV   MeV     MeV',                &
     &        '     e^2 fm^4   e^2 fm^4     e fm^2 ',            &
     &         '          e fm^2        e fm^2       ') 
  204 format (4i3,2f7.3,1f8.3,2f10.2,1f12.3,1f8.3,6f12.5)
 
!..................................................initialization
      amas = Nucl%nucleon(2)
      npro = Nucl%nucleon(1)
      r00 = 1.2d0*amas**third
      fac = 4.d0*pi/(3.d0*npro*r00**2)
      egs = GCM%EJK(0,1)   
      maxmp = GCM%NOQ
!    ............................................. E2
      call QLreduced_ME1B()
      Lif = 2 
      do ji = jmn,jmx,jdf 
         ji2   = 2*ji + 1
         kboni = GCM%kmax !min(jkmax(ji),kmx)  
         if(ji.eq.1) goto 89
            jfmax = min(ji+Lif,jmx)
         do jf = ji, jfmax, Lif 
            jf2   = 2*jf + 1
            kbonf = GCM%kmax !min(jkmax(jf),kmx)  
            fac2  = wignei(ji,2,jf,0,0,0)*sq(jf2)*iv(abs(ji-2))   
            fac22 = fac2**2   
            print 202
         do li=1,kboni 
            xei  = GCM%EJK(ji,li) - egs
         do lf=1,kbonf
            xef1 = GCM%EJK(jf,lf) - egs
!    ............................................. initialization
            be2 = zero
            be0 = zero
 
            call Calc_QLred(be0,be2,ji,jf,li,lf)
            ! be0: transform from <Jf||r^2Y_0|| Ji> to <Jf, Mf| r^2Y_00 |Ji, Mi>
            be0 = be0 / sqrt(2*ji+1.0)

            be2s  = be2**2/ji2     ! ji2 = 2*ji+1
            be22  = be2s*ji2/jf2
            beta2 = zero
!           ....................................................
            if(fac2.ne.0.d0)                                &
     &      beta2 = fac * dsqrt(be2s)/fac2   
!           .................................................
            qt    = abs(beta2)*dsqrt(16*pi/5)/fac 
            xme0  = be0     ! /sq(ji2)  is removed 2013/08/17
            rho0  = be0/(r00**2)
            rho02 = rho0**2 ! /ji2  is removed 2013/08/17
!     .............................................. print out 
          deif = xef1 - xei
          print 204,li,ji,lf,jf,xei,xef1,deif,              &
     &              be2s,be22,qt,abs(beta2),be2,xme0,rho02  
!     ..............................................
        enddo !lf          
        enddo !li   
        enddo ! jf   
   89 enddo ! ji
!
      return
      end



       subroutine Calc_QLred(q0red,q2red,ji,jf,li,lf)
       use VAPHFB_Par
       implicit none
       real*8 q0red,q2red,qred,TD1B,fact,charge
       integer tnlja,tnljb,jta,jtb,jja,Ji,Jf,Lam
       dimension TD1B(1:FSPB%tnlj_fmax,1:FSPB%tnlj_fmax,0:JJmax,0:JJmax,0:2)  !
       integer li,lf,jna,jnb,lja,ljb
       integer SPB_twoj_lj,SPB_l_lj

        do Lam =0, 2, 2
           qred = 0.d0 
        do tnlja=1,FSPB%tnlj_fmax
           jta = FSPB%t(tnlja)  ! 0, 1
           jna = FSPB%n(tnlja) ! starting from 1 
           lja = FSPB%lj(tnlja)
           jja = SPB_twoj_lj(lja)    ! doubled
           fact = 1.0/sqrt(2*Lam+1.0)
        do tnljb=1,FSPB%tnlj_fmax

           jtb = FSPB%t(tnljb)

           if(jta.ne.jtb) cycle ! does not mix n with p

           if(Input%IntType.eq.1) then
              if(jta.eq.0) charge = 0.d0 
              if(jta.eq.1) charge = 1.d0
           else
              if(jta.eq.0) charge = Input%en
              if(jta.eq.1) charge = Input%ep
           endif

           qred = qred &
                    + charge*fact  &
                    * GCMDens%TD1BSum(tnlja,tnljb,Ji,Jf,li,lf,Lam) &
                    * rME1B%T(tnlja,tnljb,Lam)

        enddo 
        enddo
          if(Lam.eq.0) q0red = qred
          if(Lam.eq.2) q2red = qred
        enddo
         
        return
        end
