!______________________________________________________________________________
      subroutine trans_by_me(jmn,jmx,jdf)

!..............................................................................
!     calculate reduced transition rates                                      .
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
  203 format (4i5,1f10.2,1f12.3,1f8.3,2f10.3)
  202 format (' i1 Ji i2 Jf   Ei    Ef    DeltaE',               &
     &        '  E2:(J_i->J_f)(J_f -> J_i) |Q(t)|',              &
     &        ' |beta(t)| <Ji||E2||Jf>   m(E0)    rho^2(E0)',/,  &
     &        '               MeV   MeV     MeV',                &
     &        '     e^2 fm^4   e^2 fm^4     e fm^2 ',            &
     &         '          e fm^2        e fm^2       ') 
  204 format (4i3,2f7.3,1f8.3,2f10.2,1f12.3,1f8.3,6f12.3)
 
!..................................................initialization
      amas = Nucl%nucleon(2)
      npro = Nucl%nucleon(1)
      r00 = 1.2d0*amas**third
      fac = 4.d0*pi/(3.d0*npro*r00**2)
      egs = GCM%EJK(0,1)   
      maxmp = GCM%NOQ
!    ............................................. E2
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
            do iqki = 1,GCM%nmaxdi(ji)  !nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = mod(iqki,maxmp)             ! mesh point only in q-space
               iki  = 2*Int((iqki-1)/maxmp)       ! k value
               if(iqi.eq.0)       iqi = maxmp
               if(iv(ji).lt.zero) iki = iki + 2   ! if jtot = 3, 5, 7, ...; k value 
            do iqkf = 1,GCM%nmaxdi(jf)                !summation over q_i and K_i 
               iqf  = mod(iqkf,maxmp)             ! mesh point only in q-space
               ikf  = 2*Int((iqkf-1)/maxmp)       ! k value
               if(iqf.eq.0)       iqf = maxmp
               if(iv(jf).lt.zero) ikf = ikf + 2   ! if jtot = 3, 5, 7, ...; k value
!    .............................................. 
               coef  =  GCM%FqJk(iqki,ji,li)*GCM%FqJk(iqkf,jf,lf) 
!    .............................................. calculations of BE0 and spectroscopic Q
               if(jf.eq.ji) then

               be2 = be2 + coef * dreal(Kernel%qpred  (iqf,iqi,ji,ikf,iki,0))
               be0 = be0 + coef * dreal(Kernel%q0p (iqf,iqi,ji,ikf,iki))
!              I found a bug here.
!              Attention: qp() is for <Jf qf || Q2 || Ji qi>
!              be2 = be2 + coef * dreal(Kernel%qp  (iqi,iqf,ji,iki,ikf))
!              be0 = be0 + coef * dreal(Kernel%q0p (iqi,iqf,ji,iki,ikf))
               endif
!    .............................................. calculations of BE2 
               if(jf.eq.ji+2)       &
     &         be2 = be2 + coef*dreal(Kernel%qpred (iqf,iqi,ji,ikf,iki,2))

!               if(jf.eq.2 .and. ji.eq.0) then
!                print *, GCM%FqJk(iqki,ji,li),GCM%FqJk(iqkf,jf,lf),dreal(Kernel%qpred (iqi,iqf,ji,iki,ikf)),be2
!               endif

               

            enddo ! iqf 
          enddo ! iqi    
            be2s  = be2**2/ji2     ! ji2 = 2*ji+1
            be22  = be2s*ji2/jf2
            beta2 = zero
!           ....................................................
            if(fac2.ne.0.d0)                                &
     &      beta2 = fac * dsqrt(be2s)/fac2   
!           .................................................
            qt    = abs(beta2)*dsqrt(16*pi/5)/fac 
            xme0  = be0     ! /sq(ji2)  is removed 2013/08/17
            rho0  = be0/r00**2
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
!______________________________________________________________________________
      subroutine trans_by_me_general(jmn,jmx,jdf)

!..............................................................................
!     calculate reduced transition rates                                      .
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
  203 format (4i5,1f10.2,1f12.3,1f8.3,2f10.3)
  202 format (' i1 Ji i2 Jf   Ei    Ef    DeltaE',               &
     &        '  E2:(J_i->J_f)(J_f -> J_i) |Q(t)|',              &
     &        ' |beta(t)| <Ji||E2||Jf>   m(E0)    rho^2(E0)',/,  &
     &        '               MeV   MeV     MeV',                &
     &        '     e^2 fm^4   e^2 fm^4     e fm^2 ',            &
     &         '          e fm^2        e fm^2       ') 
  204 format (4i3,2f7.3,1f8.3,2f10.2,1f12.3,1f8.3,6f12.3)
 
!..................................................initialization
      amas = Nucl%nucleon(2)
      npro = Nucl%nucleon(1)
      r00 = 1.2d0*amas**third
      fac = 4.d0*pi/(3.d0*npro*r00**2)
      egs = GCM%EJK(0,1)   
      maxmp = GCM%NOQ
!    ............................................. E2
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
            do iqki = 1,GCM%nmaxdi(ji)  !nmaxdi(ji)                !summation over q_i and K_i 
               iqi  = GCM%iq(iqki,ji)             ! mesh point only in q-space
               iki  = GCM%ik(iqki,ji)      ! k value
            do iqkf = 1,GCM%nmaxdi(jf)                !summation over q_i and K_i 
               iqf  = GCM%iq(iqkf,jf)             
               ikf  = GCM%ik(iqkf,jf)      ! k value
               !if(ji.eq.0 .and. jf.eq.2) print *, ji,iqi,iqf,iki,ikf,dreal(Kernel%qpred(iqf,iqi,ji,ikf,iki,2))
!    .............................................. 
               coef  =  GCM%FqJk(iqki,ji,li)*GCM%FqJk(iqkf,jf,lf) 
!    .............................................. calculations of BE0 and spectroscopic Q
               if(jf.eq.ji) then

               be2 = be2 + coef * dreal(Kernel%qpred  (iqf,iqi,ji,ikf,iki,0))
               be0 = be0 + coef * dreal(Kernel%q0p (iqf,iqi,ji,ikf,iki))


!              I found a bug here.
!              Attention: qp() is for <Jf qf || Q2 || Ji qi>
!              be2 = be2 + coef * dreal(Kernel%qp  (iqi,iqf,ji,iki,ikf))
!              be0 = be0 + coef * dreal(Kernel%q0p (iqi,iqf,ji,iki,ikf))
               endif
!    .............................................. calculations of BE2 
               if(jf.eq.ji+2)       &
     &         be2 = be2 + coef*dreal(Kernel%qpred (iqf,iqi,ji,ikf,iki,2))

!               if(jf.eq.2 .and. ji.eq.0) then
!                print *, GCM%FqJk(iqki,ji,li),GCM%FqJk(iqkf,jf,lf),dreal(Kernel%qpred (iqi,iqf,ji,iki,ikf)),be2
!               endif

               

            enddo ! iqf 
          enddo ! iqi    
            be2s  = be2**2/ji2     ! ji2 = 2*ji+1
            be22  = be2s*ji2/jf2
            beta2 = zero
!           ....................................................
            if(fac2.ne.0.d0)                                &
     &      beta2 = fac * dsqrt(be2s)/fac2   
!           .................................................
            qt    = abs(beta2)*dsqrt(16*pi/5)/fac 
            xme0  = be0     ! /sq(ji2)  is removed 2013/08/17
            rho0  = be0/r00**2
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
