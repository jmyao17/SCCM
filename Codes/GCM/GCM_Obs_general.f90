
      subroutine obser_general(jmn,jmx,jdf)
!..............................................................................
!     calculate spectroscopic quadupole moment                                . 
!..............................................................................
      USE VAPHFB_Par 
      implicit real*8(a-h,o-z)

      real*8  aa(GCM%NOQ*(JJmax*2+1),GCM%NOQ*(JJmax*2+1))
      real*8  bb(GCM%NOQ*(JJmax*2+1),GCM%NOQ*(JJmax*2+1))
      real*8  s(GCM%NOQ*(JJmax*2+1),GCM%kmax,0:JJmax)
      real*8  y(GCM%NOQ*(JJmax*2+1)),zz(GCM%NOQ*(JJmax*2+1)),w,x

      real*8, dimension(GCM%kmax,0:JJmax) ::  bet_aver, gam_aver, &
     &                 qspec,xrp,q1,q2
    
!     print *, jmn,jmx,jdf
!    ......... initialization
      q1 = zero
      q2 = zero
      maxmp = GCM%NOQ
!     .............................
      do ji = jmn,jmx,jdf
         ji2   = 2*ji + 1
         kboni = GCM%kmax ! (ji)
         if(ji.eq.1) goto 89
!     ............................................ print out (n^J(q,q))^(1/2)
         if(ji.eq.0) then
           do iq1=1,GCM%NOQ
               gg  = GCM%njkkqq(ji,iq1,iq1)
               write(23,13) iq1,sqrt(gg),GCM%FqJk(iq1,ji,1),GCM%FqJk(iq1,ji,2)
           enddo ! iq1
         endif
!     ............................................ initialization
        im   = GCM%nmaxdi(ji)
!     .................... average quadruple deformation
        do ik=1,kboni
            betav = zero
            gamav = zero
        do iqk1 = 1, im
               iq1 = GCM%iq(iqk1,ji)             ! mesh point only in q-space
               ivv   =  1
               if(abs(Const%gamma2t_mesh(iq1)-180.d0).lt.1) ivv = -1
               betav = betav + GCM%GJkq(iqk1,ji,ik)**2 * Const%beta2t_mesh(iq1)*ivv
               gamav = gamav + GCM%GJkq(iqk1,ji,ik)**2 * Const%gamma2t_mesh(iq1)
          enddo ! iqk1  
               bet_aver(ik,ji) = betav
               if(AMP%i3DAMP.eq.0) gam_aver(ik,ji) = zero
               if(AMP%i3DAMP.eq.1) gam_aver(ik,ji) = gamav
          enddo ! ik  
!       ..................................................

        do iqk1 = 1, im
        do iqk2 = 1, im
           aa(iqk1,iqk2)  = GCM%hjkkqq(ji,iqk1,iqk2) !*GCM%njkkqq(ji,iqk1,iqk2)
           bb(iqk1,iqk2)  = GCM%njkkqq(ji,iqk1,iqk2)
        enddo !  iqk1
        enddo !  iqk2
         jf   = ji
         do li= 1,kboni
            lf= li
!     ............................................ initialization
            do iqk1 = 1, im
               s(iqk1,li,ji)  = GCM%FqJk(iqk1,ji,li)
            enddo !  iqk1
            be2 = zero
            be0 = zero
            do iqki = 1,GCM%nmaxdi(ji)                !summation over q_i and K_i 
               iqi = GCM%iq(iqki,ji)          
               iki = GCM%ik(iqki,ji)      
            do iqkf = 1,GCM%nmaxdi(jf)                !summation over q_i and K_i 
               iqf = GCM%iq(iqkf,jf)       
               ikf = GCM%ik(iqkf,jf)    
!    ..............................................
               coef= GCM%FqJk(iqki,ji,li)*GCM%FqJk(iqkf,jf,lf)
               be2 = be2 + coef  * dreal(Kernel%qpred(iqf,iqi,ji,iki,ikf,0))
               be0 = be0 + coef  * dreal(Kernel%q0p (iqi,iqf,ji,iki,ikf))
            enddo ! iqf 
          enddo ! iqi    
          qspec(li,ji) = be2*dsqrt(16*pi/5)*wignei(ji,2,ji,ji,0,-ji)
          xrp(li,ji)   = dsqrt(be0)  ! /dsqrt(sq(ji2)) Aug. 17, 2013 (jmy) 
!       ........................................... accuracy of the calculation  
        ek = GCM%EJK(ji,li)
        do i=1,im
           y (i) = 0.0d0
           zz(i) = 0.0d0
           do j=1,im
              x = 0.0
!     ...........................
            do l=1,im
               x = x + (aa(j,l)-ek*bb(j,l)) * s(l,li,ji)
!                if(li.eq.1 .and. ji.eq.0) print *, s(l,li,ji)
            enddo  ! l
!     ...........................
               w    = aa(i,j)-ek*bb(i,j)
               y(i) = y(i) + w*x
!               if(li.eq.1 .and. ji.eq.0) print *, i,y(i)
               zz(i)= zz(i)+ w*s(j,li,ji)
          enddo ! j
        enddo  ! i
!     .............................
        q1(li,ji) = 0.0d0
        q2(li,ji) = 0.0d0
        do i=1,im
          q1(li,ji) = q1(li,ji) + y (i)*y (i)
          q2(li,ji) = q2(li,ji) + zz(i)*zz(i)
        enddo ! i  
!        if(li.eq.1 .and. ji.eq.0) print *, q1(li,ji), q2(li,ji)
!................................................. 
        enddo !li     
   89 enddo ! ji
!................................................. 
!
!     print out spectrum                                                      .
!                                                                             .
!     iqk = iq+k1m/2*maxmp                                                    .
!     k1m = k1-iki0, where iki0=0 (even J) or 2 (odd J)                       .
!     iq = number of points in q space                                        .
!     ik = number of points in k space                                        .
!     k  = number of eigenvalues of the norm kernel                           . 
!                                                                             .
!.............................................................................. 
!     EJK(k,J)    : energy for state J_k                                     .
!     FqJk(iqk,J,k): wave function for state J_k                              .
!..............................................................................

  13 format (i3,3f12.6)
  102 format(1x,92('-'),/)
  103 format('   Ground State energy  ',/)
  104 format('   E_g.s.= ',f10.4,' MeV')
  200 format('   excitation spectrum  ',/)
  201 format ('  J  i      E       E_ex       <beta>   <gamma>     Q_s'  &
     &  '      <N>       <Z>      rrms_p       accuracy')
  204 format (86(' '), '<H-<H>>^2')
  203 format (92(' '), '<H^2-<H>^2>^2 ')
  202 format (2i3,f12.4,3f10.4,f12.4,3f9.4,2f9.4)
!..................................................initialization
      print 102
      amas = Nucl%nucleon(2)
      r00 = 1.2d0*amas**third
      egs = GCM%EJK(0,1)
      print 103
      write(*,104) egs
      print 102
      print 200
      print 201
      print 203
      print 204
      print 102

      maxmp = GCM%NOQ 
!    ...........................................
      do ji   = jmn,jmx,jdf
         ji2  = 2*ji + 1
         if(ji.eq.1) goto 99
            kbon  = GCM%kmax !min(jkmax(ji),kmx)
         do li    = 1,kbon
            xei   = GCM%EJK(ji,li) - egs
!    ............................................. initialization 
            qtot  = zero
            xn    = zero
            xp    = zero
!           ftest = zero
            do iqki = 1,GCM%nmaxdi(ji)                !summation over q_i and K_i 

               iqi = GCM%iq(iqki,ji)             
               iki = GCM%ik(iqki,ji) 

            do iqkf = 1,GCM%nmaxdi(ji)                !summation over q_i and K_i 

               iqf = GCM%iq(iqkf,ji)             
               ikf = GCM%ik(iqkf,ji) 

               coef = GCM%FqJk(iqki,ji,li)*GCM%FqJk(iqkf,ji,li)
!    .............................................. calculations of <N> & <Z>  
              xn = xn + coef * dreal(Kernel%nn(iqi,iqf,ji,iki,ikf))  &
     &                       * GCM%njkkqq(ji,iqkf,iqki)
              xp = xp + coef * dreal(Kernel%zz(iqi,iqf,ji,iki,ikf)) &
     &                       * GCM%njkkqq(ji,iqkf,iqki)

!            ftest = ftest + coef*GCM%njkkqq(ji,iqkf,iqki) 
            enddo ! iqkf 
          enddo ! iqki
!          write(*,*) 'ftest=',ftest 
          xrp(li,ji) = xrp(li,ji)/dsqrt(xp)          ! normalized to proton number
          write(*,202) ji, li, GCM%EJK(ji,li), xei, bet_aver(li,ji),  &
     &              gam_aver(li,ji),qspec(li,ji),xn,xp,xrp(li,ji), &
     &              q1(li,ji),q2(li,ji)
        enddo !li    
   99 enddo ! ji
!
      print 102
      return
      end

