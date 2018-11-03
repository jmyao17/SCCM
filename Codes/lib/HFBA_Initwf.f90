       subroutine Init_hfbwf()
       use VAPHFB_Par
       implicit none

       integer ii,jj


      if(.NOT. ALLOCATED(HFB%U0))  ALLOCATE(HFB%U0(1:HO%NLEV,1:HO%NLEV))
      if(.NOT. ALLOCATED(HFB%V0))  ALLOCATE(HFB%V0(1:HO%NLEV,1:HO%NLEV))


       print *, ' -> Initial WF '
        if(Input%inwf.eq.0) then
         do ii=1,HO%NLEV
          do jj=1,HO%NLEV
           read(10,*) HFB%U0(ii,jj),HFB%V0(ii,jj)
          end do
         end do
         call checkwf(HFB%U0,HFB%V0,HO%NLEV)
         close(10)
        elseif (Input%inwf.eq.1) then
         print*,'    Read the init w.f. from save file ...'
         call initwf2(HFB%U0,HFB%V0,HO%NLEV)

        else if(Input%inwf.eq.2) then
         print*,'    another arbitrary init w.f.'
         call initwf2(HFB%U0,HFB%V0,HO%NLEV)
         print*,'    check init w.f.'
         call checkwf(HFB%U0,HFB%V0,HO%NLEV)
        else
         print*,'    Invalid inwf !!!';stop
        end if

        return
        end
!      ....................................
       subroutine initwf(U0,V0,jshort,NDIM)
!      ....................................
       USE VAPHFB_PAR
       USE IFPORT
       implicit real*8(a-h,o-z)

        DIMENSION U0(NDIM,NDIM),V0(NDIM,NDIM)
        DIMENSION AAA(NDIM,NDIM)
	DIMENSION BBB(NDIM,NDIM)
	DIMENSION AUX1(NDIM,NDIM)
	DIMENSION AUX2(NDIM,NDIM)
	DIMENSION AUX3(NDIM,NDIM)
	DIMENSION AUX4(NDIM,NDIM)
	DIMENSION WEVA(NDIM)
	DIMENSION WORKA(3*NDIM-1)
	DIMENSION WEVB(NDIM)
	DIMENSION WORKB(3*NDIM-1)
	DIMENSION jshort(3)
	
	
!POINTERS
	do ii=1,ndim
	 do jj=1,ndim
           U0(ii,jj)=zero
           V0(ii,jj)=zero
	  end do
	end do
	
	ipos=0
	do i=1,3
	j=jshort(i)
	mmax=(j+1)/2
	do m=1,mmax
	ipos1=ipos+m
	ipos2=ipos+j+2-m
	ipos3=ipos1+NDIM/2
	ipos4=ipos2+NDIM/2

	 U0(ipos1,ipos1)= dsqrt(1.d0*rand())
	 U0(ipos2,ipos2)= U0(ipos1,ipos1)

	 U0(ipos3,ipos3)=  dsqrt(1.d0*rand())
	 U0(ipos4,ipos4)= U0(ipos3,ipos3)

	 V0(ipos1,ipos2)= dsqrt(1.d0-(U0(ipos1,ipos1))**2)
	 V0(ipos2,ipos1)=-V0(ipos1,ipos2) 

	 V0(ipos3,ipos4)= dsqrt(1.d0-(U0(ipos3,ipos3))**2)  
	 V0(ipos4,ipos3)=-V0(ipos3,ipos4)


	end do
	ipos=ipos+j+1
	end do
       
	end subroutine



        subroutine initwf2(U0,V0,NDIM)
        use IFPORT

         implicit real*8(a-h,o-z)
         Parameter (zero=0.d0)
         Parameter (one=1.d0)
         DIMENSION U0(NDIM,NDIM),V0(NDIM,NDIM)
         DIMENSION AAA(NDIM,NDIM)
         DIMENSION BBB(NDIM,NDIM)
         DIMENSION AUX1(NDIM,NDIM)
         DIMENSION AUX2(NDIM,NDIM)
         DIMENSION AUX3(NDIM,NDIM)
         DIMENSION AUX4(NDIM,NDIM)
         DIMENSION WEVA(NDIM)
         DIMENSION WORKA(3*NDIM-1)
         DIMENSION WEVB(NDIM)
         DIMENSION WORKB(3*NDIM-1)

         print *, ' -> initial wf2'
!POINTERS
         ip  = 0
         ipb = ip+NDIM/4
         in  = ipb+NDIM/4
         inb = in+NDIM/4        
 
         do ii=1,ndim
          do jj=1,ndim
             U0(ii,jj)=zero
             V0(ii,jj)=zero
          end do
        end do

         print *, ' generating random numbers for (U,V)'
        do ii=1,NDIM/4
           U0(ii+ip,ii+ip )  = dsqrt(1.d0*rand())
           U0(ii+ipb,ii+ipb )= U0(ii+ip,ii+ip )
           U0(ii+in,ii+in )  = dsqrt(1.d0*rand())
           U0(ii+inb,ii+inb )= U0(ii+in,ii+in )
           V0(ii+ip,ii+ipb ) = dsqrt(1.d0-(U0(ii+ip,ii+ip ))**2)
           V0(ii+ipb,ii+ip ) =-V0(ii+ip,ii+ipb ) 
           V0(ii+in,ii+inb ) = dsqrt(1.d0-(U0(ii+in,ii+in ))**2)
  
           V0(ii+inb,ii+in )= -V0(ii+in,ii+inb )
        end do

        do ii=1,NDIM
          do jj=ii,NDIM
             AAA(jj,ii)=rand()
             BBB(jj,ii)=rand()
          end do
        end do

         do ii=1,NDIM
         do jj=1,NDIM
            AAA(ii,jj)=AAA(jj,ii)
            BBB(ii,jj)=BBB(jj,ii)
          end do
          end do



        call dsyev('v','L',NDIM,AAA,NDIM,WEVA,WORKA,3*NDIM-1,INFO)
        call dsyev('v','L',NDIM,BBB,NDIM,WEVB,WORKB,3*NDIM-1,INFO)
         
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,AAA,NDIM,U0,NDIM,zero,AUX1,NDIM)	
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,AUX1,NDIM,BBB,NDIM,zero,AUX2,NDIM)

        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,AAA,NDIM,V0,NDIM,zero,AUX3,NDIM)	
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,AUX3,NDIM,BBB,NDIM,zero,AUX4,NDIM)

        do ii=1,NDIM
          do jj=1,NDIM
             U0(jj,ii)= AUX2(jj,ii)
             V0(jj,ii)= AUX4(jj,ii)
          end do
        end do
        end subroutine


        subroutine checkwf(U0,V0,NDIM)
        USE VAPHFB_PAR
	implicit real*8(a-h,o-z)

	DIMENSION U0(NDIM,NDIM),V0(NDIM,NDIM)
	DIMENSION AUX1(NDIM,NDIM)
	DIMENSION AUX2(NDIM,NDIM)
	DIMENSION AUX3(NDIM,NDIM)
	DIMENSION AUX4(NDIM,NDIM)	
	DIMENSION AUX5(NDIM,NDIM)
	DIMENSION AUX6(NDIM,NDIM)	
	DIMENSION BOGO1(NDIM,NDIM)
	DIMENSION BOGO2(NDIM,NDIM)
	
	
        open(2,file='initwf.dat',status='unknown')
!	do ii=1,NDIM
!	do jj=1,NDIM
!	  write(2,'(2i5,2f12.8)') ii,jj,U0(ii,jj),V0(ii,jj)
!	end do
!	end do
!	write(2,*),'***********************'
	
	call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,U0,NDIM,U0,NDIM,0.d0,AUX1,NDIM)	
        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,V0,NDIM,V0,NDIM,0.d0,AUX2,NDIM)
  

        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,U0,NDIM,V0,NDIM,0.d0,AUX3,NDIM)
        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,V0,NDIM,U0,NDIM,0.d0,AUX4,NDIM)
        call DGEMM ('n','t',NDIM,NDIM,NDIM,-1.d0,U0,NDIM,V0,NDIM,0.d0,AUX6,NDIM)

        do ii=1,NDIM
	  do jj=1,NDIM
	    BOGO2(ii,jj)=AUX3(ii,jj)+AUX4(ii,jj)
	    BOGO1(ii,jj)=AUX1(ii,jj)+AUX2(ii,jj)
	  end do
	end do
	
        write(2,*),' Checking the relation: U^+U + V^+V = I'
        do ii=1,NDIM
        do jj=1,NDIM

        if(abs(BOGO1(ii,jj)).gt. 0.1) write(2,'(2i5,2f12.8)') ii,jj,BOGO1(ii,jj)
        end do
        end do
        write(2,*),'***********************'
        write(2,*),' '
        return

        write(2,*),'UV++V*Ut'
	do ii=1,NDIM
	  write(2,1000) (BOGO2(ii,jj), jj=1,NDIM)
	end do
	write(2,*),'***********************'
	write(2,*),' '


	write(2,*),'Density Matrix'
	do ii=1,NDIM
	  write(2,1000) (AUX2(ii,jj), jj=1,NDIM)
	end do          
	write(2,*),'***********************'
	write(2,*),' '


	write(2,*),'Pairing Tensor 10'
	do ii=1,NDIM
	  write(2,1000) (AUX3(ii,jj), jj=1,NDIM)
	end do  


1000    format(10000F20.8)

        END SUBROUTINE
