!______________________________________________________________________________
          subroutine  norm_overlap(npz2,znorm)
!______________________________________________________________________________
!    determination of overlap with the Pfafiian method
!    with the formula by Bertsch and Robledo, PRL108, 042505 (2012).
!
!    1. construct the skew-symmetric matrix M (4n * 4n)
!      where 2n is the dimension of the matrix V and U.
!
!          (  V(a)^T U^(a) ,  V(a)^T R^T V^{(b)*} )
!      M = (   ------------   ------------------  )
!          (               ,  U^{(b)+} V{(b)^*}   )
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
      integer i,j,ir,jr,npz2,npz4
      complex*16 znorm
!      parameter (ntest=10)
      character*1 UPLO, MTHD
      complex*16 Ma(npz2,npz2),Mb(npz2,npz2),Mc(npz2,npz2)
      complex*16 cva,cvb
      complex*16 Zmax(npz2*2,npz2*2)
      complex*16 ipiv(npz2*2,npz2*2)
      complex*16 ztmp

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
      write(*,*) ' Matrix a'
        call zmaxa(npz2,Ma(1,1),cva)
      write(*,*) ' Matrix b'
        call zmaxb(npz2,Mb(1,1),cvb)
      write(*,*) ' Matrix c'
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
      write(*,*) ' Pfaffian'
        call ZPfaffianF(Zmax(1,1),npz4,npz4,ipiv(1,1),ztmp)
!------------------------------------------------------
      write(*,*) 'cnorm=',cva*cvb,'ztmp1=',ztmp  
      znorm = iv(npz2/2)*ztmp/(cva*cvb)
!     write(*,*) '<R> by Pf=',znorm
!----END NORM
      return
      END
!........................................................................
      subroutine zmaxa(npz2,Ma,cna)
!........................................................................
!     construct the M^(a)
!     in a canoical basis
!     1 <-> a
!     2 <-> b
!     <a | b>
!----------------------------------------------------------------
      use VAPHFB_Par
      implicit real*8 (a-h,o-z)
      complex*16 Ma,zmatr
      dimension Ma(npz2,npz2),zmatr(npz2,npz2)

!    ....................................................................
         do i=1,npz2    ! dimension of full model space, including (k, kbar)
            do j=1,npz2
                zmatr(i,j)=zero
            enddo ! j
          enddo ! i
!     ............... normalization factor
          cna = one
          if(iv(npz2).lt.0) stop 'npz2 is not an even number !!'
          do k = 1,npz2, 2
!     ................................................... partially
!     occupied s.p.state, n and p states are not separated
             if(HFB%vv1(k).le.PNP%eps1) stop 'error in truncation'  ! -> q2
                vk  = dsqrt(HFB%vv1(k))
                uk  = dsqrt(1.d0-HFB%vv1(k))
!     ..........................................matrix elements of f
!     (representation in s.p.states)
                zmatr(k,k+1)  = -vk*uk
                zmatr(k+1,k)  =  vk*uk
                cna  =  cna * vk
!          write(*,*) k,vk,uk
  3       enddo !k
!     ........................................................ 
       do i=1,npz2
       do j=1,npz2
           Ma(j,i) = zmatr(j,i)
       enddo ! j
       enddo ! i
!       do i=1,npz2
!       do j=i,npz2
!         if(abs(Ma(j,i)).gt.1.d-7) 
!       write(*,'(2i3,4f15.8)') i,j,Ma(i,j),Ma(j,i)
!       enddo ! i
!       enddo ! i
      return
      end
!==============================================================================
!........................................................................
      subroutine zmaxb(npz2,Mb,cnb)
!........................................................................
!     construct the M^(b)
!     in a canoical basis
!     1 <-> a
!     2 <-> b
!     <a | b>
!----------------------------------------------------------------
      use VAPHFB_Par
      implicit real*8 (a-h,o-z)
      complex*16 Mb,zmatr
      dimension Mb(npz2,npz2),zmatr(npz2,npz2)

!    ....................................................................
         do i=1,npz2
            do j=1,npz2
                zmatr(i,j)=zero
            enddo ! j
          enddo ! i
          cnb = one
!     ...........................

          if(iv(npz2).lt.0) stop 'npz2 is not an even number !!'
          do k = 1, npz2, 2
!     ................................................... partially
!     occupied state
             if(HFB%vv2(k).le.PNP%eps2) stop 'error in truncation'  ! -> q2
                vk  = dsqrt(HFB%vv2(k))
                uk  = dsqrt(1.d0-HFB%vv2(k))
!     ..........................................matrix elements of f
!     (representation in s.p.states)
                zmatr(k,k+1)  =  vk*uk
                zmatr(k+1,k)  = -vk*uk
                cnb  =  cnb * vk
  3       enddo !k
!     ........................................................ 
       do i=1,npz2
       do j=1,npz2
           Mb(j,i) = zmatr(j,i)
       enddo ! j
       enddo ! i

      return
      end


!-------------------------------------------------------------------------------
      subroutine zmaxab(npz2,Mc)
!........................................................................
!     construct M^(c)
!     in a canoical basis
!........................................................................
      use vaphfb_par
      implicit real*8 (a-h,o-z)
      complex*16 Rmatrix,rot,rotive,Rive,on
      complex*16 zphi_n,zphi_p,zeiphi,zeiphi_n,zeiphi_p,zei2phi_n,zei2phi_p
      complex*16 zmatr,zmatrT
      complex*16 Mc,Rmat(npz2,npz2),eiphi
      dimension zmatr(npz2,npz2),Mc(npz2,npz2),zmatrT(npz2,npz2)
      complex*16 ZAUX1(HO%NLEV,HO%NLEV),ZAUX2(HO%NLEV,HO%NLEV) 
      complex*16 vmatr,vmatrT
      dimension vmatr(npz2,npz2),vmatrT(npz2,npz2)
      complex*16 ipiv(npz2,npz2),zdetR
!    ....................................................................
      zphi_n    = zimag*PNP%phi(1)
      zeiphi_n  = cdexp(zphi_n)
!      zei2phi_n = cdexp(2*zphi_n)
      zphi_p    = zimag*PNP%phi(2)
      zeiphi_p  = cdexp(zphi_p)
!      zei2phi_p = cdexp(2*zphi_p)

!    ...........................
          do i=1,npz2
          do j=1,npz2
             zmatr(i,j) = zero
             zmatrT(i,j) = zero
             vmatr(i,j) = zero
             vmatrT(i,j) = zero
          enddo ! j
          enddo ! i

!     ...........................
          if(iv(npz2).lt.0) stop 'npz2 is not an even number !!'
          do k = 1, npz2,2
!     ................................................... partially occupied state
             if(HFB%vv1(k).le.PNP%eps1) stop 'error in trucation' 
!     ....... check the occupation prob. for time-reveral states
             if(abs(HFB%vv1(k)-HFB%vv1(k+1)) .gt.1.d-4) then 
               write(*,*) HFB%vv1(k),HFB%vv1(k+1)
               stop 'error in vv'
             endif
                vk = dsqrt(HFB%vv1(k))
                uk = dsqrt(1.d0-HFB%vv1(k))
!     ..........................................matrix elements of V^{(a)T}
!     (representation in s.p.states)
                vmatrT(k,k+1)   = -vk
                vmatrT(k+1,k)   =  vk
  2       enddo !k



          if(iv(npz2).lt.0) stop 'npz2 is not an even number !!'
          do k = 1,npz2,2

             if(HFB%it2(k).gt.0) zeiphi = zeiphi_n 
             if(HFB%it2(k).lt.0) zeiphi = zeiphi_p 
!     ................................................... partially
!     occupied state
             if(HFB%vv2(k).le.PNP%eps2) stop 'error in trucation' 
                vk = dsqrt(HFB%vv2(k))
                uk = dsqrt(1.d0-HFB%vv2(k))
!     ..........................................matrix elements of V^{(b)}
!     (representation in s.p.states)
                vmatr(k,k+1)  =  zeiphi*vk
                vmatr(k+1,k)  = -zeiphi*vk
          enddo !k

!     ........................................................
!     computing R(Omega)_jk
!     ........................................................
           do j=1,npz2
           do k=1,npz2
              ZAUX2(j,k)=zzero
           do m1=1,HO%NLEV
           do m2=1,HO%NLEV
              ZAUX2(j,k)=ZAUX2(j,k) + HFB%RO_1(m1,j)*AMP%ZROT_m1m2(m1,m2)*HFB%RO_2(m2,k)
           enddo ! m1
           enddo ! m2 
!           if(abs(ZAUX2(j,k)).gt.1.d-7) write(*,'(2i4,2f10.5)') j,k,ZAUX2(j,k) 
!           if(abs(cmatr(j,k)).gt.1.d-7) write(*,'(2i4,2f10.5)') j,k,cmatr(j,k) 
           enddo ! k
           enddo ! j
 
!     ........................................................
                do j=1,npz2
                do k=1,npz2
                   zmatr(j,k) = zero
                do l=1,npz2
                   zmatr(j,k) = zmatr(j,k)  + ZAUX2(l,j) * vmatr(l,k)
                enddo ! l
                enddo ! k
                enddo ! j


                do i=1,npz2
                do j=1,npz2
                   Mc(i,j) = zero
                do l=1,npz2
                   Mc(i,j) = Mc(i,j)  + vmatrT(i,l)*zmatr(l,j)
                enddo ! l
                enddo ! j
                enddo ! i
  10    format(40f20.18)
!     ........................................................
      return
      end

