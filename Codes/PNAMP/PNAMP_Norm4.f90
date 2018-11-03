!    ...................................................
        subroutine  norm_overlap(npz2,znorm)
!    ...................................................
!    Neergard-Wuest formulas
!    NPA 402, 311 (1983)
!    ...................................................
      use vaphfb_par
      implicit none
      integer INFO,npz22,npz2,npz4,i,j,ir,jr
      complex*16 znorm,ztmp
!     ..............................................................      
      complex*16 zphi_n,zphi_p,zei2phi,zeiphi_i,zeiphi_j,zeiphi_n,zeiphi_p
      Integer N,m1,m2,k,l
      complex*16 ZZTZ(npz2,npz2),Zb_bar(npz2,npz2),Zb(npz2,npz2),Za(npz2,npz2) 
      real*8     test,vl,ul,vk,uk,vkl,ukl
 
      complex*16 ZWORK2,ZWORK1,ZWR,ZVL,ZVR
      real*8    cM(npz2,npz2),WR(npz2),em(npz2),Rmat(npz2,npz2),cna,cnb,WORK1
      Dimension ZWR(npz2)
      Dimension ZVL(npz2,npz2)
      Dimension ZVR(npz2,npz2)
      Dimension ZWORK1(2*npz2)
      Dimension ZWORK2(npz2)
      Dimension WORK1(2*npz2)
!     ........................................................
!     computing D = R(Omega)_jk
!     ........................................................
           write(*,*) ' computing D ..' 
           do j=1,npz2
           do k=1,npz2
              Rmat(j,k) = zzero
           do m1=1,HO%NLEV
           do m2=1,HO%NLEV
              Rmat(j,k)=Rmat(j,k) + HFB%RO_1(m1,j)*AMP%ZROT_m1m2(m1,m2)*HFB%RO_2(m2,k)
           enddo ! m1
           enddo ! m2 
           enddo ! k
           enddo ! j

           Za = zzero
           Zb = zzero
           cna = one
           cnb = one 

!    ... attention to the ordering of energy levels, time-reversal states
!     Za: 
!     computing C' = (B'A'^-1)* 
!                  = ( 0       v1/u1  0      0       ...  )
!                    (-v1/u1    0     0      0       ...  )
!                    (  0       0     0      v2/u2   ...  )
!                    (  0       0    -v2/u2  0       0    )
!                    ( ....    ....                       )
!    ...................................
!     Zb_bar = D Zb D^T: 
!     computing C  = (BA^-1)* 
!                  = ( 0       v1/u1  0      0       ...  )
!                    (-v1/u2    0     0      0       ...  )
!                  D (  0       0     0      v2/u2   ...  ) D^T
!                    (  0       0    -v2/u2  0       0    )
!                    ( ....    ....                       )
!    ...................................
          do i=1,npz2,2
             ir = i+1
             vk = dsqrt(HFB%vv1(i))
             uk = dsqrt(1.d0-HFB%vv1(i))
             vl = dsqrt(HFB%vv2(i))
             ul = dsqrt(1.d0-HFB%vv2(i))
             if(abs(uk).lt.1.d-12) stop 'uk is close to zero'
             if(abs(ul).lt.1.d-12) stop 'ul is close to zero'
             Za(i,ir) =  vk/uk
             Za(ir,i) = -vk/uk

             Zb(i,ir) =  vl/ul
             Zb(ir,i) = -vl/ul
 
             cna = cna * uk
             cnb = cnb * ul
           enddo ! i

       do i=1,npz2
       do j=1,npz2
          Zb_bar(i,j) = zzero
       do k=1,npz2
       do l=1,npz2
          Zb_bar(i,j) = Zb_bar(i,j) + Rmat(i,k)*Zb(k,l)*Rmat(j,l) 
        enddo ! i
        enddo ! i
!        write(*,'(2i3,4f12.6)') i,j,Za(i,j),Zb_bar(i,j)
        enddo ! i
        enddo ! i


       do i=1,npz2
       do j=1,npz2
          cM(i,j) = zero
       do k=1,npz2
          cM(i,j) = cM(i,j) + real(Za(k,i)*Zb_bar(k,j)) 
        enddo ! i
!        if(i.eq.1 .and. j.eq.4) write(*,'(2i3,f12.6)') i,j,cM(i,j)
        enddo ! i
        enddo ! i

      call sdiag(npz2,npz2,cM,WR,cM,em,1)
!      do i=1,npz2
!          write(*,*) i,WR(i)
!      enddo

      znorm=zone
      do i=1,npz2,2  ! eigenvalues are pair-wise
       znorm = znorm * (one+WR(i))
      end do
      znorm  = cna*cnb*znorm   ! two-fold degenerate
!     ............................................................
      write(*,*) 'znorm=',znorm,cna*cnb
      return
      end
!........................................................................

 
