!       ...................................
        subroutine Isospin_t3()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

      if(.NOT. ALLOCATED(t3me))  ALLOCATE(t3me(1:HO%NLEV,1:HO%NLEV))

        do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             t3me(jj,ii)=zero
             if((ii.eq.jj).AND.(ii.le.(HO%NLEV/2))) then
                t3me(jj,ii)=-1.d0
             else if((ii.eq.jj).AND.(ii.gt.(HO%NLEV/2))) then
                t3me(jj,ii)=1.d0
             end if
          end do
       end do
       end subroutine

!       ...................................
        subroutine ZNme()
!       ...................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
 
      if(.NOT. ALLOCATED(AZme))  ALLOCATE(AZme(1:HO%NLEV*HO%NLEV))
      if(.NOT. ALLOCATED(ANme))  ALLOCATE(ANme(1:HO%NLEV*HO%NLEV))

       ll=0
        do ii=1,HO%NLEV
          do jj=1,HO%NLEV
             ll=ll+1
             AZme(ll)=zero
             ANme(ll)=zero
             if((ii.eq.jj).AND.(ii.le.(HO%NLEV/2))) then
                AZme(ll)=one
             else if((ii.eq.jj).AND.(ii.gt.(HO%NLEV/2))) then
                ANme(ll)=one
             end if
          end do
       end do
       end subroutine



       subroutine JZ_JZ2_ME(AJX,AJY,AJZ,AJX2,AJY2,AJZ2,NLEV)
     
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)

        DIMENSION AJX(NLEV,NLEV)
        DIMENSION AJX2(NLEV,NLEV)
        DIMENSION AJY(NLEV,NLEV)
        DIMENSION AJY2(NLEV,NLEV)
        DIMENSION AJZ(NLEV,NLEV)
        DIMENSION AJZ2(NLEV,NLEV)
        b2=HO%b_osc**2
!initializing
         do ii=1,HO%NLEV
         do jj=1,HO%NLEV
            AJX(jj,ii)=zero
            AJX2(jj,ii)=zero
            AJY(jj,ii)=zero
            AJY2(jj,ii)=zero
            AJZ(jj,ii)=zero
            AJZ2(jj,ii)=zero
        end do
        end do

        do ia=1,NLEV
         ina=tnljm%n(ia) ! nlindex(ia)
         ila =tnljm%l(ia) ! lang(ia)
         ija =tnljm%twoj(ia) ! jang(ia)
         ima =tnljm%twom(ia) ! mjang(ia)
         ita =tnljm%t(ia)    ! mtisos(ia)

         rja=ija/2.d0
         rma=ima/2.d0

         do ib=1,NLEV
         inb=tnljm%n(ib) ! nlindex(ia)
         ilb =tnljm%l(ib) ! lang(ia)
         ijb =tnljm%twoj(ib) ! jang(ia)
         imb =tnljm%twom(ib) ! mjang(ia)
         itb =tnljm%t(ib)    ! mtisos(ia)

         rjb=ijb/2.d0
         rmb=imb/2.d0

 
         fp1=rjb*(rjb+1.d0)-rmb*(rmb+1.d0)
         fm1=rjb*(rjb+1.d0)-rmb*(rmb-1.d0)

         fp1p2=rjb*(rjb+1.d0)-(rmb+1.d0)*(rmb+2.d0)
         fm1m2=rjb*(rjb+1.d0)-(rmb-1.d0)*(rmb-2.d0)
         f02=  rjb*(rjb+1.d0)-rmb**2

         if ((ina.eq.inb).AND.(ila.eq.ilb)            &
     &   .AND.(ija.eq.ijb).AND.(ita.eq.itb)) then
         
         delta_ma_mbp1=0.d0
         delta_ma_mbm1=0.d0
         delta_ma_mb  =0.d0
         delta_ma_mbp2=0.d0
         delta_ma_mbm2=0.d0

         if ((ima.eq.(imb+2))) delta_ma_mbp1=1.d0
         if ((ima.eq.(imb-2))) delta_ma_mbm1=1.d0
         if ((ima.eq.(imb+4))) delta_ma_mbp2=1.d0
         if ((ima.eq.(imb-4))) delta_ma_mbm2=1.d0
         if ((ima.eq.imb)) delta_ma_mb=1.d0
 
         AJX(ia,ib)=0.5d0*(dsqrt(fp1)*delta_ma_mbp1+               &
     &                    dsqrt(fm1)*delta_ma_mbm1)  
         AJY(ia,ib)=0.5d0*(dsqrt(fp1)*delta_ma_mbp1-               &
     &                     dsqrt(fm1)*delta_ma_mbm1)
     
         AJZ(ia,ib)=rmb*delta_ma_mb
     
         AJX2(ia,ib)=0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+       &
     &                       dsqrt(fm1*fm1m2)*delta_ma_mbm2+         &
     &                                2*f02*delta_ma_mb)

         AJY2(ia,ib)=0.25d0*(dsqrt(fp1*fp1p2)*delta_ma_mbp2+       &
     &                      dsqrt(fm1*fm1m2)*delta_ma_mbm2-       &
     &                                2*f02*delta_ma_mb)
         
             AJZ2(ia,ib)=(rmb*delta_ma_mb)**2
           end if
          end do
        end do
        end subroutine
!      ..........................................
      subroutine Pairs_TMT_JMJ(IT,MMT,JJ,MMJ,P_TMT_JMJ_me,NLEV)
!     ....................................
!       iso-scalar pairing and iso-vector pairing
!       with the pairing amplitude:
!       P_TMT_JMJ =: 1/sqrt(2) * sum_l (2*l+1)^1/2 [c^+_l c^+l]^LST,
!       where L=0, S=1, T=0 for iso-scalar pairing 
!             L=0, S=0, T=1 for iso-vector pairing 
!     ....................................
!     In JT-scheme:
!     iso-vector pairing (S=0,T=1) => (J=0,T=1)
!     P_TMT_JMJ = 1/2 * sum_j (2*j+1)^1/2 [c^+_j c^+j]^(JT)_(0mu),
!               = 1/2 * sum_(jm) sum_(t,t') <1/2 t, 1/2 t | T=1, mu>
!                     * C^+_jmt C^+_jmt'                
!     ....................................

        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        DIMENSION P_TMT_JMJ_me(NLEV,NLEV)
        REAL(DP) :: CG
        do ia=1,NLEV
          jla=tnljm%l(ia) 
          jja=tnljm%twoj(ia) !jang(ia)
          jma=tnljm%twom(ia) !mjang(ia)
          jta=tnljm%t(ia)    !mtisos(ia)

          fact=dsqrt(jja+1.d0)/2.d0   ! 1/2 * (2*j+1)^1/2
         do ib=1,NLEV
!       ............. initialization
           P_TMT_JMJ_me(ia,ib)= 0.d0

          jlb=tnljm%l(ib) 
          jjb=tnljm%twoj(ib) !jang(ia)
          jmb=tnljm%twom(ib) !mjang(ia)
          jtb=tnljm%t(ib)    !mtisos(ia): -1 (p) and +1 (n)

          delta_ab=0.d0
          cb1 = 0.d0
          cb3 = 0.d0
!      .......... iso-vector pairing case (T=1, S=0?)
          if(iT.eq.1 .and. JJ.eq.0) then
             if(jja.eq.jjb) delta_ab=1.d0
!             call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)
!             call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)

             cb1 = CG(jja,jma,jjb,jmb,2*JJ,2*MMJ) 
             cb3 = CG(1,jta,1,jtb,2*IT,2*MMT) 
             P_TMT_JMJ_me(ia,ib)=fact*cb1*cb3*delta_ab

!      .......... iso-scalar pairing case
          elseif(iT.eq.0 .and. JJ.eq.1) then
             if(jla.eq.jlb) delta_ab=1.d0         ! delta_(la,lb)
!             call CJJ(jja,jjb,2*JJ,jma,jmb,2*MMJ,cb1)  
!             call CJJ(1,1,2*IT,jta,jtb,2*MMT,cb3)      
             cb1 = CG(jja,jma,jjb,jmb,2*JJ,2*MMJ) 
             cb3 = CG(1,jta,1,jtb,2*IT,2*MMT) 
             P_TMT_JMJ_me(ia,ib) = 2*delta_ab*cb1*cb3  &   ! multiplied by a factor 2 by hand
     &         *Wigner_6j(jjb,jja,2*JJ,1,1,2*jla)    &
     &         *iv(jla+(jja+3)/2)*dsqrt((jja+1.d0)*(jjb+1.d0))/dsqrt(2.d0)
          else
             stop ' Error in the definition of pairing !'
          endif
         end do
        end do
                
        end subroutine


!     .....................................
        subroutine Pairs_TLS(IT,MMT,IS,MMS,P_TMT_SMS_me,NLEV)
!     ....................................
!       Iso-scalar pairing amplitude:
!       P_ST = 1/sqrt(2) * sum_l (2*l+1)^1/2 [c^+_l c^+l]^LST,
!       where L=0, S=1, T=0 for iso-scalar pairing 
!             L=0, S=0, T=1 for iso-vector pairing 
!     ....................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        integer IS
        Real(DP) :: CG
        DIMENSION P_TMT_SMS_me(NLEV,NLEV)

        do ia=1,NLEV
          jla=tnljm%l(ia) 
          jja=tnljm%twoj(ia) !jang(ia)
          jma=tnljm%twom(ia) !mjang(ia)
          jta=tnljm%t(ia)    !mtisos(ia)

!         fact=dsqrt(jja*1.d0+0.5d0)
!          fact=dsqrt(jja/2.d0+0.5d0)
         do ib=1,NLEV

          jlb=tnljm%l(ib) 
          jjb=tnljm%twoj(ib) !jang(ia)
          jmb=tnljm%twom(ib) !mjang(ia)
          jtb=tnljm%t(ib)    !mtisos(ia)

          delta_ab=0.d0
          if(jla.eq.jlb) delta_ab=1.d0

          P_TMT_SMS_me(ia,ib) = zero
 
          call CJJ(jja,jjb,IS*2,jma,jmb,2*MMS,cbj)       ! isospin part
          call CJJ(1,1,2*IT,jta,jtb,2*MMT,cbt)       ! isospin part

          P_TMT_SMS_me(ib,ia) = P_TMT_SMS_me(ia,ib) &
     &        + delta_ab*cbj*cbt*Wigner_6j(jjb,jja,2*IS,1,1,2*jla) &
     &         *iv(jla+(jja+3)/2)*dsqrt((jja+1.d0)*(jjb+1.d0))/dsqrt(2.d0)
!      ..............
        enddo
        enddo
        return
        end subroutine


        SUBROUTINE Qpair_0(QME,akapa10,akapa01,Q_0,NDIM)
        
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
!        Parameter (zero=0.d0)
!        Parameter (one=1.d0)
        
        DIMENSION QME(NDIM,NDIM)
        DIMENSION akapa10(NDIM,NDIM),akapa01(NDIM,NDIM)
        DIMENSION AUX1(NDIM,NDIM)
        DIMENSION AUX2(NDIM,NDIM)

        
        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,    &
     &  QME,NDIM,akapa10,NDIM,zero,AUX1,NDIM)   

        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,    &
     &  QME,NDIM,akapa01,NDIM,zero,AUX2,NDIM)   
        
        
        Q_0=zero
        do ii=1,NDIM
         Q_0=Q_0+(AUX1(ii,ii)+AUX2(ii,ii))/2.d0
        end do
        
        END SUBROUTINE
