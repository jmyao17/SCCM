!______________________________________________________________________________
          subroutine  norm_overlap(npz2,znorm)
!______________________________________________________________________________
!    determination of overlap with the Pfafiian method
!    with the formula by Bertsch and Robledo, PRL108, 042505 (2012).
!
!    1. construct the skew-symmetric matrix M (4n * 4n)
!      where 2n is the dimension of the matrix V and U.
!
!          (  V(a)^T U^(a) ,  V(a)^T V^{(b)*}     )
!      M = (   ------------   ------------------ )
!          (               ,  U^{(b)+} V{(b)^+}  )
!
!          (      M^(a) ,           M^(c)        )
!        = (   ------------   ------------------ )
!          (               ,        M^(b)        )
!
!    2. calculate the Pfaffian of M
!
!    3. calculate the overlap <qa|qb> = (-1)^n * Pf(M)/ (N^(a)*N^(b))
!       where normalization factor N^(a) = Prod_i v^(a)_i.
!
!     calculate M matrice in norm overlap: <q_a|R|q_b>
!                                 NOTICE : <q_1|R|q_2>  
!    4. return the norm overlap with phase
!----------------------------------------------------------------

!----------------------------------------------------------------
      use vaphfb_par
      implicit none
      integer npz2,npz4
      complex*16 znorm
!      parameter (ntest=10)
      character*1 UPLO, MTHD
      complex*16 Ma(npz2,npz2),Mb(npz2,npz2),Mc(npz2,npz2)
      complex*16 cva,cvb,pcva,pcvb
      complex*16 Zmax(npz4,npz4)
      complex*16 ipiv(npz4,npz4)
      complex*16 ztmp

!      complex*16 AV(npz4,npz4),EL(npz4),WORK(npz4*3)

!      dimension IWORK(npz4),RWORK(npz4)
!      complex*16 WORK2(npz4)
!       dimension IWORK_t(ntest),RWORK_t(ntest)
!       complex*16 aa(ntest,ntest),WORK2_t(ntest),
!     &            iaa(ntest,ntest)

      data UPLO/'U'/,MTHD/'P'/
!     ............................................................. initialization
       npz4 = npz2*2
       do i=1,npz4
       do j=1,npz4
           Zmax(j,i) = zero
       enddo ! j
       enddo ! i
!     ..............................................................      
!     cna: prod_k u^a_k;  cnb: prod_l u^b_l
!     ..............................................................      
        call zmaxa(npz2,Ma(1,1),cva)
        call zmaxb(npz2,Mb(1,1),cvb)
        call zmaxab(npz2,Mc(1,1))
!     .......................... construction of matrix M
       do i=1,npz2
          ir = i + npz2
       do j=1,npz2
          jr = j + npz2
          Zmax(i,j)    =  Ma(i,j)
          Zmax(i,jr)   =  Mc(i,j)
          Zmax(jr,i)   = -Mc(i,j)
          Zmax(ir,jr)  =  Mb(i,j)
       enddo ! j
       enddo ! i 
!------------------------------------------------------
!     to calculate the Pfaffian of a skew-symmetric matrix  
!------------------------------------------------------
!       TWO SUBROUTINES for Pfaffian of a complex matrix
!    The second one is (10 times) faster than the first one
!------------------------------------------------------

        call ZPfaffianF(Zmax(1,1),npz4,npz4,ipiv(1,1),ztmp)
!------------------------------------------------------
!        call ZSKPFA(UPLO,MTHD,npz4,Zmax,npz4,ztmp,
!     $              IWORK,WORK2,npz4, RWORK, INFO)
!------------------------------------------------------
!      write(*,*) 'cnorm=',cna*cnb,'ztmp1=',ztmp  
!------------------------------------------------------
       on = iv(npz2/2)*ztmp/(cva*cvb)
!------------------------------------------------------

!      write(*,*) 'it=',it,'cnorm=',cna*cnb,'ztmp=',ztmp  
     write(*,*) '<R> by Pf=',on
!----END NORM
      return
      END
!-------------------------------------------------------------------------------
      subroutine zmaxab(npz2,Mc)
!........................................................................
!     construct M^(c)
!     in a canoical basis
!........................................................................
      implicit real*8 (a-h,o-z)
      complex*16 Rmatrix,rot,rotive,Rive,on
      complex*16 zphi_n,zphi_p,zei2phi_n,zei2phi_p
      complex*16 zmatr,zmatrT
      complex*16 zmatr,Mc,Rmat(npz2,npz2),eiphi
      dimension zmatr(npz2,npz2),Mc(npz2,npz2),zmatr(npz2,npz2),zmatrT(npz2,npz2)
!    ....................................................................
      zphi_n    = cmplx(0,PNP%phi(1))
      zeiphi_n  = cdexp(zphi_n)
      zei2phi_n = cdexp(2*zphi_n)
      zphi_p    = cmplx(0,PNP%phi(2))
      zeiphi_p  = cdexp(zphi_p)
      zei2phi_p = cdexp(2*zphi_p)

          do i=1,npz2
          do j=1,npz2
             zmatr(i,j) = zero
             zmatrT(i,j) = zero
          enddo ! j
          enddo ! i
!     ...........................
          do k = 1, npz2 
!     ................................................... partially occupied state
             if(HFB%vv1(k).le.PNP%eps1) goto 2
                vk = dsqrt(HFB%vv1(k))
                uk = dsqrt(1.d0-HFB%vv1(k))
                ik  = HFB%kocc1(k)

!                ikr = ik+npz
!                ik2 = 2*ik
!     ..........................................matrix elements of V^{(a)T}
!     (representation in s.p.states)
                cmatrT(ik2-1,ik2)   = -vk
                cmatrT(ik2,ik2-1)   =  vk
  2       enddo !k
          k1 = ka2(ib,it) + 1
          k2 = ka2(ib,it) + kd2(ib,it)
          do k = k1,k2
!     ................................................... partially
!     occupied state
             if(HFB%vv2(k).le.HFB%eps2) goto 3
                vk = dsqrt(HFB%vv2(k))
                uk = dsqrt(1.d0-HFB%vv2(k))
                ik  = HFB%kocc2(k)
                ikr = ik+npz(ib,it)
                ik2 = 2*ik
!     ..........................................matrix elements of V^{(b)}
!     (representation in s.p.states)
                cmatr(ik2-1,ik2)  =  eiphi*vk
                cmatr(ik2,ik2-1)  =  -eiphi*vk
  3       enddo !k

         npz1 = npz(ib,it)
!           write(97,*) '----- R----' 
         do i = 1,npz1
            ir= i + npz1
            i2= i*2
         do j = 1,npz1
            jr= j + npz1
            j2= j*2
            Rmat(i2-1,j2-1) = Rmatrix(i,j,ib,it)
            Rmat(i2-1,j2  ) = Rmatrix(i,jr,ib,it)
            Rmat(i2  ,j2-1) = Rmatrix(ir,j,ib,it)
            Rmat(i2  ,j2  ) = Rmatrix(ir,jr,ib,it)
         enddo !i
!           write(97,10) (abs(Rmat(i,j)),j=1,npz2)
         enddo !i
!           write(97,*) '----- Mc----' 
!     ........................................................
                do j=1,npz2
                do k=1,npz2
                   zmatr(j,k) = zero
                do l=1,npz2
                   zmatr(j,k) = zmatr(j,k)    + Rmat(j,l) * cmatr(l,k)
                enddo ! l
                enddo ! k
                enddo ! j
                do i=1,npz2
                do j=1,npz2
                   Mc(i,j) = zero
                do l=1,npz2
                   Mc(i,j) = Mc(i,j)  + cmatrT(i,l)*zmatr(l,j)
                enddo ! l
                enddo ! j
!                  write(97,10) (abs(Mc(i,j)),j=1,npz2)
                enddo ! i
  10    format(40f20.18)
!     ........................................................
      return
      end

!........................................................................
      subroutine zmaxa(npz2,Ma,cna)
!........................................................................
!     construct the M^(a)
!     in a canoical basis
!     1 <-> a
!     2 <-> b
!     <a | b>
!----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 Ma,cmatr
      dimension Ma(npz2,npz2),cmatr(npz2,npz2)

!    ....................................................................
         do i=1,npz2
            do j=1,npz2
                cmatr(i,j)=zero
            enddo ! j
          enddo ! i
          cna = one
!     ...........................
          k1 = ka1(ib,it) + 1
          k2 = ka1(ib,it) + kd1(ib,it)
          do k = k1,k2
!     ................................................... partially
!     occupied state
             if(HFB%vv1(k).le.PNP%eps1) goto 3   ! -> q2
                vk  = dsqrt(HFB%vv1(k))
                uk  = dsqrt(1.d0-HFB%vv1(k))
                ik  = HFB%kocc1(k)
                ikr = ik+npz(k,it)
!     ..........................................matrix elements of f
!     (representation in s.p.states)
                ik2 = ik*2
                cmatr(ik2-1,ik2)  = -vk*uk
                cmatr(ik2,ik2-1)  =  vk*uk
                cna  =  cna * vk
  3       enddo !k
!     ........................................................ 
       do i=1,npz2
       do j=1,npz2
           Ma(j,i) = cmatr(j,i)
       enddo ! j
       enddo ! i

      return
      end
c==============================================================================
!........................................................................
      subroutine zmaxb(npz2,Mb,cnb)
!........................................................................
!     construct the M^(b)
!     in a canoical basis
!     1 <-> a
!     2 <-> b
!     <a | b>
!----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 Mb,cmatr
      dimension Mb(npz2,npz2),cmatr(npz2,npz2)

!    ....................................................................
         do i=1,npz2
            do j=1,npz2
                cmatr(i,j)=zero
            enddo ! j
          enddo ! i
          cnb = one
!     ...........................

          k1 = ka2(ib,it) + 1
          k2 = ka2(ib,it) + kd2(ib,it)
          do k = k1,k2
!     ................................................... partially
!     occupied state
             if(vocc2(k,it).le.epsocc2(ib,it)) goto 3   ! -> q2
                vk  = dsqrt(vv2(k,it))
                uk  = dsqrt(1.d0-vv2(k,it))
                ik  = HFB%kocc2(k,it)
                ikr = ik+npz(ib,it)
!     ..........................................matrix elements of f
!     (representation in s.p.states)
                ik2 = ik*2
                cmatr(ik2-1,ik2)  =  vk*uk
                cmatr(ik2,ik2-1)  = -vk*uk
                cnb  =  cnb * vk
  3       enddo !k
!     ........................................................ 
       do i=1,npz2
       do j=1,npz2
           Mb(j,i) = cmatr(j,i)
       enddo ! j
       enddo ! i

      return
      end

