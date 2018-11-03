      subroutine Model_Space_Generator(emax,nljmax)

      implicit none
      integer NN,emax,n,l,l2,j2,nlj,nljmax
      integer emax1,emax2,it

      emax1 =mod(emax/10,10) + 48
      emax2 =mod(emax,10) + 48


!#....... Model Space..... order should be consistent with int.
!    nlj    t      n      l    twoj
!Neutron:  (newly added)

      open(7,file='emax'//char(emax1)//char(emax2)//'.val',status='unknown')

      write(7,*) "....... Model Space..... order should be consistent with int."
      write(7,*) "    nlj    t      n      l    twoj "
! ...........................................
      do it=0,1
         if(it.eq.0) write(7,*) "Neutron"
         if(it.eq.1) write(7,*) "Proton"
         nlj=0
      do NN=0,eMax
         do L=0,NN
            n=(NN-L)/2
            L2 = L*2
            if(2*n+L.ne.NN) cycle
            do j2=L2+1,max(1,L2-1),-2
               nlj = nlj+1
               write(7,'(5i8)') nlj,it,n,l,j2
            enddo 
         enddo 
      enddo 
      enddo
      nljmax= nlj*2
      return
      end
      
