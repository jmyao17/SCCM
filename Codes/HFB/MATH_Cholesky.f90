! .................................................................................................
!   +       t *
! LL = 1 + Z Z
!
! the Cholesky decomposition or Cholesky factorization is a decomposition of a Hermitian, 
! positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose
! .................................................................................................
	SUBROUTINE Lmatrix(Z,GLOB_LOC,NDIM)
	USE VAPHFB_PAR
	implicit real*8 (a-h,o-z)
	
	DIMENSION Z(NDIM,NDIM)	
	DIMENSION GLOB_LOC(NDIM,NDIM)
	DIMENSION AUX1(NDIM,NDIM)
	DIMENSION AUX2(NDIM,NDIM)
	DIMENSION AUX3(NDIM+1,NDIM)
	

	call DGEMM ('t','n',NDIM,NDIM,NDIM,one,             &
     &	Z,NDIM,Z,NDIM,zero,AUX1,NDIM)

 
     	do ii=1,NDIM
	  do jj=ii,NDIM
	    if (ii.ne.jj) then
	    AUX2(ii,jj)=AUX1(ii,jj)
	    AUX2(jj,ii)=AUX2(ii,jj)
	    else if (ii.eq.jj) then
	    AUX2(ii,jj)=AUX1(ii,jj)+1.d0
	    AUX2(jj,ii)=AUX2(ii,jj)
	    end if
	  end do
	end do
 
	call cholesky_format1(AUX2,AUX3,NDIM) 
	
	call DPBTRF( 'L', NDIM, NDIM,AUX3,NDIM+1, INFO )
	if (info.ne.0) then 
	 print*,'Cholesky fails.'
	 print*,'Try another seed or lower eta1,eta2'
	 stop
	end if
	
	call cholesky_format2(AUX3,GLOB_LOC,NDIM)
1000	format(1000F10.4)
	END SUBROUTINE Lmatrix


	SUBROUTINE cholesky_format1(A,AB,NDIM)
	
	implicit real*8 (a-h,o-z)
	DIMENSION A(NDIM,NDIM)
	DIMENSION AB(NDIM+1,NDIM)

        do ii=1,NDIM
	  do jj=1,NDIM
	    kk=1+ii-jj
	    if(jj.le.ii) then
	      AB(kk,jj)=A(ii,jj)
	    end if
	  end do
        end do


	END SUBROUTINE

	SUBROUTINE cholesky_format2(AB,A,NDIM)
	
	implicit real*8 (a-h,o-z)
	
	PARAMETER (one=1.d0)
	PARAMETER (zero=0.d0)
	DIMENSION A(NDIM,NDIM)
	DIMENSION AB(NDIM+1,NDIM)

	
	 do ii=1,NDIM
	   do jj=1,NDIM
	     A(jj,ii)=zero
	   end do
	 end do
	
        do ii=1,NDIM
	  do jj=1,NDIM+1-ii
         
	   ll = (ii-1)+jj
	   mm = jj
          
	   A(ll,mm)=AB(ii,jj)
          
	  end do
 
        end do


	 
	
	END SUBROUTINE
	


	SUBROUTINE Lmatrix_VAP(Z,GLOB_LOC,NDIM)
        USE VAPHFB_PAR
	implicit real*8 (a-h,o-z)
	
	COMPLEX*16 IM,Izero,Ione
	PARAMETER (IM=(0.d0,1.d0))
	PARAMETER (Izero=(0.d0,0.d0))
	PARAMETER (Ione=(1.d0,0.d0))
	
	COMPLEX*16 Z(NDIM,NDIM)	
	COMPLEX*16 GLOB_LOC(NDIM,NDIM)
	COMPLEX*16 AUX1(NDIM,NDIM)
	COMPLEX*16 AUX2(NDIM,NDIM)
	COMPLEX*16 AUX3(NDIM+1,NDIM)
	

	call ZGEMM ('c','n',NDIM,NDIM,NDIM,Ione,   &
     &	Z,NDIM,Z,NDIM,Izero,AUX1,NDIM)

 
     	do ii=1,NDIM
	  do jj=1,NDIM
	    if (ii.ne.jj) then
	    AUX2(jj,ii)=dreal(AUX1(jj,ii))-IM*dimag(AUX1(jj,ii))
	    else if (ii.eq.jj) then
	    AUX2(jj,ii)=dreal(AUX1(jj,ii))-IM*dimag(AUX1(jj,ii))+ one
	    end if
	  end do
	end do
	call cholesky_format1_VAP(AUX2,AUX3,NDIM) 
	
	call ZPBTRF( 'L', NDIM, NDIM,AUX3,NDIM+1, INFO )
	if (info.ne.0) then 
	 print*,'Cholesky fails.'
	 print*,'Try another seed or lower eta1,eta2'
	 stop
	end if
	
	call cholesky_format2_VAP(AUX3,GLOB_LOC,NDIM)
1000	format(1000F10.4)
	END SUBROUTINE 


	
	SUBROUTINE cholesky_format1_VAP(A,AB,NDIM)
	
	implicit real*8 (a-h,o-z)
	COMPLEX*16 A(NDIM,NDIM)
	COMPLEX*16 AB(NDIM+1,NDIM)

        do ii=1,NDIM
	  do jj=1,NDIM
	    kk=1+ii-jj
	    if(jj.le.ii) then
	      AB(kk,jj)=A(ii,jj)
	    end if
	  end do
        end do


	END SUBROUTINE

	SUBROUTINE cholesky_format2_VAP(AB,A,NDIM)
	
	implicit real*8 (a-h,o-z)
	COMPLEX*16 Izero
	PARAMETER (Izero=(0.d0,0.d0))
	COMPLEX*16 A(NDIM,NDIM)
	COMPLEX*16 AB(NDIM+1,NDIM)

	
	 do ii=1,NDIM
	   do jj=ii,NDIM
	     A(ii,jj)=Izero
	     A(jj,ii)=Izero
	   end do
	 end do
	
        do ii=1,NDIM
	  do jj=1,NDIM+1-ii
         
	   ll = (ii-1)+jj
	   mm = jj
          
	   A(ll,mm)=AB(ii,jj)
          
	  end do
 
        end do

	END SUBROUTINE
