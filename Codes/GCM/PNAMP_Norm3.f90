!    ...................................................
        subroutine  norm_overlap(npz2,znorm)
!    ...................................................
      use vaphfb_par
      implicit none
      integer npz22,npz2,npz4,i,j,ir,jr
      complex*16 znorm,ztmp
!     ..............................................................      
      complex*16 zphi_n,zphi_p,zei2phi,zeiphi_i,zeiphi_j,zeiphi_n,zeiphi_p
      Integer m1,m2,k,l
      complex*16 DetR,Rmat(npz2,npz2),CDET,Dmat(npz2,npz2),zrot(1:npz2,1:npz2) 
      complex*16 Rive(npz2,npz2)
      real*8     test,vl,ul,vk,uk,vkl,ukl
!    ....................................................................
      zphi_n    = zimag*PNP%phi(1)
      zeiphi_n = cdexp(zphi_n)
      zphi_p    = zimag*PNP%phi(2)
      zeiphi_p = cdexp(zphi_p)
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
!       do m1=1,HO%NLEV !npz2
!       do m2=1,HO%NLEV !npz2
!          write(*,'(2i3,4f12.6)') m1,m2, AMP%ZROT_m1m2(m1,m2),AMP%ZROT_m1m2(m2,m1) !Rmat(i,j),Rmat(j,i)
!        enddo ! i
!        enddo ! i
!        stop

      DETR = zone
      call rmatrixive(npz2,Rmat,DETR,Rive)
      write(*,*) 'DetR=',DETR
!     ................................................... partially occupied state

            do i=1,npz2,2
               ir = i+1
               if(ir.gt.npz2) stop 'error'
            do j=1,npz2,2
               jr = j+1
               if(jr.gt.npz2) stop 'error'
             if(HFB%it1(i).gt.0) zeiphi_i = zeiphi_n
             if(HFB%it1(i).lt.0) zeiphi_i = zeiphi_p
             if(HFB%it2(j).gt.0) zeiphi_j = zeiphi_n
             if(HFB%it2(j).lt.0) zeiphi_j = zeiphi_p
             zei2phi = zeiphi_i*zeiphi_j
!    ... attention to the ordering of energy levels, time-reversal states
             vk = dsqrt(HFB%vv1(i))
             vl = dsqrt(HFB%vv2(j))
             uk = dsqrt(1.d0-HFB%vv1(i))
             ul = dsqrt(1.d0-HFB%vv2(j))
            vkl = vk*vl
            ukl = uk*ul
       zrot(i,j)   = DCONJG(Rmat(i,j))*vkl*zei2phi  + Rive(i,j)*ukl
       zrot(i,jr)  = DCONJG(Rmat(i,jr))*vkl*zei2phi + Rive(i,jr)*ukl
       zrot(ir,j)  = DCONJG(Rmat(ir,j))*vkl*zei2phi + Rive(ir,j)*ukl
       zrot(ir,jr) = DCONJG(Rmat(ir,jr))*vkl*zei2phi+ Rive(ir,jr)*ukl

           enddo ! j
           enddo ! i

!       do i=1,npz2
!       do j=1,npz2
!          write(*,'(2i3,4f12.6)') i,j, Rmat(i,j),Rmat(j,i)
!        enddo ! i
!        enddo ! i
!     .................................. calculate the inverse of N: rotive and its determinant 
       CDET = zone 
       call normive(npz2,zrot,CDET)
!       call rmatrixive(npz2,zrot,CDET,Rive) 
      write(*,*) 'DetD=',CDET
!     ............................................................
        znorm  = cdsqrt(CDET * DETR)
        write(*,*) 'znorm=',znorm

      return
      end
!........................................................................

!______________________________________________________________________________
      subroutine rmatrixive(npz2,Rmat,detr,Rive)
!=======================================================================
!     calculate the inverse of R^T matrix and its determinant       
!-----------------------------------------------------------------------
      use vaphfb_par
      implicit real*8 (a-h,o-z)
      complex*16 a,detr,x,d,rot,rotive,Rmat(npz2,npz2),Rive
      dimension  a(npz2*npz2),x(npz2*npz2),Rive(npz2,npz2)
!----------------------------
       ndm=npz2
            do il=1,ndm
            do ih=1,ndm
!   ........................ R^T
             a(ih+(il-1)*ndm)=Rmat(il,ih)
             x(ih+(il-1)*ndm)=zzero
            enddo !il
            x(il+(il-1)*ndm)=zone !CMPLX(1.d0,0.d0)
            enddo !ih

            call clingd(ndm,ndm,ndm,ndm,a,x,d,ifl)

!---- Rive is the inverse of R^T
!            do il=1,ndm
!            do ih=1,ndm
!               Rive(ih,il)=x(ih+(il-1)*ndm)
!            enddo !il 
!            enddo !ih

            detr=detr*d

      return
      END

!______________________________________________________________________________
      subroutine normive(ndm,zrot,cdet)
!=======================================================================
!     calculate the inverse of D matrix and its determinant       
!-----------------------------------------------------------------------
      use vaphfb_par
      implicit real*8 (a-h,o-z)
      complex*16 za,cdet,zx,d,zrot
      dimension  za(1:ndm*ndm),zx(1:ndm*ndm),zrot(1:ndm,1:ndm)

        write(*,*) ' normive: ndm=',ndm
            do il=1,ndm
            do ih=1,ndm
               za(ih+(il-1)*ndm) = zrot(ih,il)
               zx(ih+(il-1)*ndm) = zzero
            enddo !il
               zx(il+(il-1)*ndm) = zone ! CMPLX(1.d0,0.d0)
            enddo !ih

!            do il=1,ndm
!            do ih=1,ndm
!               write(*,'(2i4,2f12.6)') ih,il,za(ih+(il-1)*ndm)
!            enddo !ih
!            enddo !ih
        call clingd(ndm,ndm,ndm,ndm,za,zx,d,ifl)
        write(*,*) ' normive: clingd completed'

         cdet=cdet*d

!         do il=1,ndm
!            do ih=1,ndm
!                AMP%rotive(ih,il)=x(ih+(il-1)*ndm)
!           enddo !il 
!         enddo !ih
      return
!-END norminve
      END 
 
