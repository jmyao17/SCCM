!    ...................................................
        subroutine  norm_overlap(npz2,znorm)
!    ...................................................
!    determination of overlap with the Pfafiian method
!    introduced by L. M. Robledo, PRC79, 021302 (R) (2009)
!
!    1. construct the Matrix M
!          (  M^(a) ,  -1      )
!      M = (  ------   ------- )
!          (  +1    ,  -M^(b)* )
!
!    2. calculate the Pfaffian of M
!
!    3. calculate the overlap <q2|q1> = S(N) * Pf(M)
!       where S(N)=(-1)^(N(N+1)/2).
!
!     calculate M matrice in norm overlap: <q_a|R|q_b>
!                                 NOTICE : <q_1|R|q_2>  
!    4. return the norm overlap with phase
!----------------------------------------------------------------
      use vaphfb_par
      implicit none
      integer npz22,npz2,npz4,i,j,ir,jr
      real*8 sn
      complex*16 znorm,ztmp,cna,cnb,M1(npz2,npz2),M2(npz2,npz2),Zmax(npz2*2,npz2*2)
      integer ipiv(npz2*2,npz2*2)
!     ..............................................................      
       npz4 = npz2*2
       npz22 = npz2*(npz2+1)/2
       if(npz22.le.igfv) then
         sn   = iv(npz22)
        else
!          stop ' igfv is too small !'
         sn   = (-1)**npz22
       endif
!     ............................................................. initialization
       do i=1,npz4
       do j=1,npz4
           Zmax(j,i) = zzero
       enddo ! j
       enddo ! i
!     ..............................................................      
!     cna: prod_k u^a_k;  cnb: prod_l u^b_l
!     ..............................................................      
        call zrmax(npz2,M2(1,1),cnb)
        call zlmax(npz2,M1(1,1),cna)
!     .......................... construction of matrix M
       do i=1,npz2
          ir = i + npz2
          Zmax(i,ir)   = -zone
          Zmax(ir,i)   =  zone
       do j=1,npz2
          jr = j + npz2
          Zmax(j,i)   = M2(j,i)
          Zmax(jr,ir) = -dconjg(M1(j,i))
       enddo ! j
       enddo ! i
!------------------------------------------------------
!     to calculate the Pfaffian of a skew-symmetric matrix  
!------------------------------------------------------
!       TWO SUBROUTINES for Pfaffian of a complex matrix
!    The second one is (10 times) faster than the first one
!------------------------------------------------------
       call ZPfaffianF(Zmax(1,1),npz4,npz4,ipiv(1,1),ztmp)
!        call ZSKPFA(UPLO,MTHD,npz4,Zmax,npz4,ztmp,
!     $              IWORK,WORK2,npz4, RWORK, INFO)
!------------------------------------------------------
       znorm = sn*ztmp*cna*cnb
!      znorm = sn*ztmp !*cna*cnb
      write(*,*) 'cnorm=',cna*cnb,'ztmp=',ztmp  
!     write(*,*) '<R> by Pf=',on
!----END NORM
      return
      END

!    ....................................................................
      subroutine zrmax(npz2,Mb,cnb)
      use vaphfb_par
      implicit none
      complex*16 zphi_n,zphi_p,zei2phi,zei2phi_n,zei2phi_p
      integer m1,m2,i,j,npz2,k,l
      complex*16 cmatr(npz2,npz2),Mb(npz2,npz2),cnb,zmatr(npz2,npz2) 
      real*8    vk,uk
      complex*16 test,Rmat(npz2,npz2)
!    ....................................................................
      zphi_n    = zimag*PNP%phi(1)
      zei2phi_n = cdexp(2*zphi_n)
      zphi_p    = zimag*PNP%phi(2)
      zei2phi_p = cdexp(2*zphi_p)

          do i=1,npz2
          do j=1,npz2
             cmatr(i,j) = zzero
          enddo ! j
          enddo ! i
          cnb = zone
!     ...........................
          do k = 1,npz2,2
!     ................................................... partially occupied state
             if(HFB%it2(k).gt.0) zei2phi = zei2phi_n
             if(HFB%it2(k).lt.0) zei2phi = zei2phi_p

!      .................................
             if(HFB%vv2(k).le.PNP%eps2) cycle 
             if(abs(1.d0-HFB%vv2(k)).le.1.d-12) then
                write(*,*) 'v2(',k,')=',HFB%vv2(k)
                write(*,*) 'u2 is very close to 0'
                stop
             endif
!      .................................
                vk = dsqrt(HFB%vv2(k))
                uk = dsqrt(1.d0-HFB%vv2(k))
!     ..........................................matrix elements of f (representation in s.p.states)
                cmatr(k,k+1)  = zei2phi*vk/uk
                cmatr(k+1,k)  =-zei2phi*vk/uk
                cnb = cnb * uk
  3       enddo !k
!     ........................................................
!     ........................................................
!     computing R(Omega)_jk
!     ........................................................
           do j=1,npz2
           do k=1,npz2
              Rmat(j,k) = zzero
           do m1=1,HO%NLEV
           do m2=1,HO%NLEV
              Rmat(j,k)=Rmat(j,k) + HFB%RO_1(m1,j)*AMP%ZROT_m1m2(m1,m2)*HFB%RO_2(m2,k)
           enddo ! m1
           enddo ! m2 
!           if(abs(Rmat(j,k)).gt.1.d-5) write(*,'(2i3,2f12.6)') j,k,Rmat(j,k)
           enddo ! k
           enddo ! j
!     ........................................................
!                do i=1,npz2
!                do k=1,npz2
!                   test = zero   
!                   do j=1,npz2
!
!                   test = test + Rmat(i,j)*Rmat(k,j) 
!
!                enddo ! j
!                   if(abs(test).gt.0.01) write(*,*) i,k,test
!                enddo ! k
!                enddo ! i
!          stop
!     ........................................................
                do j=1,npz2
                do k=1,npz2
                   zmatr(j,k) = zzero
                do l=1,npz2
                   zmatr(j,k) = zmatr(j,k) + Rmat(j,l) * cmatr(l,k)
                enddo ! l
                enddo ! k
                enddo ! j
!     ........................................................
                do i=1,npz2
                do j=1,npz2
                   Mb(i,j) = zzero
                do l=1,npz2
                   Mb(i,j) = Mb(i,j) + zmatr(i,l) * Rmat(j,l)
                enddo ! l
                enddo ! j
                enddo ! i
!     ........................................................
      return
      end
!........................................................................
      subroutine zlmax(npz2,Ma,cna)
!........................................................................
!     construct the left (q1) Thouless matrix Z   
!     in a canoical basis
!     1 <-> a
!     2 <-> b
!     <b | a>
!----------------------------------------------------------------
      use vaphfb_par
!      implicit real*8 (a-h,o-z)
      integer npz2,i,j,k
      real*8  vk,uk
      complex*16 cna,Ma,cmatr
      dimension Ma(npz2,npz2),cmatr(npz2,npz2)
!    ....................................................................
         do i=1,npz2
            do j=1,npz2
                cmatr(i,j)=zzero
            enddo ! j
          enddo ! i
          cna = zone
!     ...........................
          do k =1,npz2,2
!     ................................................... partially occupied state
             if(HFB%vv1(k).le.PNP%eps1) cycle   ! -> q2 
             if(abs(1.d0-HFB%vv1(k)).le.1.d-12) then
                write(*,*) 'v2(',k,')=',HFB%vv1(k)
                write(*,*) 'u2 is very close to 0'
                stop
             endif
                vk  = dsqrt(HFB%vv1(k))
                uk  = dsqrt(1.d0-HFB%vv1(k))
!     ..........................................matrix elements of f (representation in s.p.states)
                cmatr(k,k+1)  = vk/uk
                cmatr(k+1,k)  =-vk/uk
                cna = cna * uk
  3       enddo !k
!     ........................................................
       do i=1,npz2
       do j=1,npz2
           Ma(j,i) = cmatr(j,i)
       enddo ! j
       enddo ! i

      return
      end
      
