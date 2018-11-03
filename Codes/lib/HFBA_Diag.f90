      subroutine HFB_Diag(hh,Delta,N) 
      use VAPHFB_Par
      implicit none
      integer N,nhfb 
      real*8  hh(N,N),Delta(N,N)
      real*8  fg(N*2*N*2),HM(N*2*N*2)
      real*8  ez(N*2),e(N*2),d(N*2)
      real*8  vvequi(N),eeequi(N),deequi(N)
      real*8  ala,sn
      integer i,j,k,klp

!   ........
!   hh does not include the constraint terms
!   ........
      nhfb = 2*N
      do i=1, N
      do j=1, N
         HM(  i+  (j-1)*nhfb)     = hh(i,j)   
         HM(N+i+(N+j-1)*nhfb)     =-hh(i,j)   
         HM(N+i+(j-1)*nhfb)       = Delta(i,j)
         HM(N+j+(i-1)*nhfb)       = Delta(i,j)
       enddo 
         if(tnljm%t(i).eq.1)  ala= HFB%EFermi(0) ! n
         if(tnljm%t(i).eq.-1) ala= HFB%EFermi(1) ! p
        HM(i+(j-1)*nhfb) = HM(i+(j-1)*nhfb) - ala
        HM(N+i+(N+j-1)*nhfb) = HM(N+i+(N+j-1)*nhfb) + ala
       enddo 

     call sdiag(N,N,HM,e,HM,ez,+1)
     ! print out eigenvalue
      do i=1,nhfb
         write(110,*) e(i)
      enddo  

      return
      end 
