

        subroutine HFB_20(U,V,ham,delta10,H20,NDIM)
!      ....................................................
!
!      20    +       *    +   t   *    +       *   +   *    *
!     H   = U * h * V  - V * h * U  + U * d * U - V * d  * V
!
!      ....................................................
        
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        
        DIMENSION U(NDIM,NDIM)
        DIMENSION V(NDIM,NDIM)
        DIMENSION ham(NDIM,NDIM)
        DIMENSION H20(NDIM,NDIM)
        DIMENSION H20AUX(NDIM,NDIM)
        DIMENSION delta10(NDIM,NDIM)
        DIMENSION AUX1(NDIM,NDIM)
        DIMENSION AUX2(NDIM,NDIM)
        DIMENSION AUX3(NDIM,NDIM)
        DIMENSION AUX4(NDIM,NDIM)
        DIMENSION AUX5(NDIM,NDIM)
        DIMENSION AUX6(NDIM,NDIM)
        DIMENSION AUX7(NDIM,NDIM)
        DIMENSION AUX8(NDIM,NDIM)
        

        
        call DGEMM ('t','n',NDIM,NDIM,NDIM,one,     &
     &	U,NDIM,ham,NDIM,zero,AUX1,NDIM)	
            call DGEMM ('n','n',NDIM,NDIM,NDIM,one, &
     &	AUX1,NDIM,V,NDIM,zero,AUX2,NDIM)	
         
            call DGEMM ('t','t',NDIM,NDIM,NDIM,one, &
     &	V,NDIM,ham,NDIM,zero,AUX3,NDIM)	
            call DGEMM ('n','n',NDIM,NDIM,NDIM,one, &
     &	AUX3,NDIM,U,NDIM,zero,AUX4,NDIM)	
         
            call DGEMM ('t','n',NDIM,NDIM,NDIM,one, &
     &	U,NDIM,delta10,NDIM,zero,AUX5,NDIM)	
            call DGEMM ('n','n',NDIM,NDIM,NDIM,one, &
     &	AUX5,NDIM,U,NDIM,zero,AUX6,NDIM)	
         
            call DGEMM ('t','n',NDIM,NDIM,NDIM,one, &
     &	V,NDIM,delta10,NDIM,zero,AUX7,NDIM)	
            call DGEMM ('n','n',NDIM,NDIM,NDIM,one, &
     &	AUX7,NDIM,V,NDIM,zero,AUX8,NDIM)

        do ii=1,NDIM
          do jj=1,NDIM
            H20AUX(jj,ii)= AUX2(jj,ii)-AUX4(jj,ii)+  &
     &	    AUX6(jj,ii)-AUX8(jj,ii) 
          end do
        end do
        do ii=1,NDIM
          do jj=1,NDIM
            H20(jj,ii)= 0.5*(H20AUX(jj,ii)  -H20AUX(ii,jj))
          end do
        end do
        
        return
        END SUBROUTINE  
               
         


        subroutine HFB_GRAD_PHI(ZU,ZV,Zham,Zdelta10,zdelta01,ZH20,NDIM)
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)
        DIMENSION ZU(NDIM,NDIM)
	DIMENSION ZV(NDIM,NDIM)
	DIMENSION Zham(NDIM,NDIM)
	DIMENSION ZH20(NDIM,NDIM)
	DIMENSION Zdelta10(NDIM,NDIM)
	DIMENSION Zdelta01(NDIM,NDIM)
	DIMENSION ZAUX1(NDIM,NDIM)
	DIMENSION ZAUX2(NDIM,NDIM)
	DIMENSION ZAUX3(NDIM,NDIM)
	DIMENSION ZAUX4(NDIM,NDIM)
	DIMENSION ZAUX5(NDIM,NDIM)
	DIMENSION ZAUX6(NDIM,NDIM)
	DIMENSION ZAUX7(NDIM,NDIM)
	DIMENSION ZAUX8(NDIM,NDIM)
	

	!Ut h V
	call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,  &
     &	ZU,NDIM,Zham,NDIM,zzero,ZAUX1,NDIM)	
     	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &	ZAUX1,NDIM,ZV,NDIM,zzero,ZAUX2,NDIM)	
     
	!Vt ht U
     	call ZGEMM ('t','t',NDIM,NDIM,NDIM,zone,  &
     &	ZV,NDIM,Zham,NDIM,zzero,ZAUX3,NDIM)	
     	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &	ZAUX3,NDIM,ZU,NDIM,zzero,ZAUX4,NDIM)	
     
	!Ut Delta10 U
     	call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone,  &
     &	ZU,NDIM,Zdelta10,NDIM,zzero,ZAUX5,NDIM)	
     	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,  &
     &	ZAUX5,NDIM,ZU,NDIM,zzero,ZAUX6,NDIM)	
     
	!Vt Delta01 V
     	call ZGEMM ('t','n',NDIM,NDIM,NDIM,zone, &
     &	ZV,NDIM,Zdelta01,NDIM,zzero,ZAUX7,NDIM)	
     	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone, &
     &	ZAUX7,NDIM,ZV,NDIM,zzero,ZAUX8,NDIM)

	do ii=1,NDIM
	  do jj=1,NDIM
	    ZH20(jj,ii)= ZAUX2(jj,ii)-ZAUX4(jj,ii)+  &
     &	    ZAUX6(jj,ii)-ZAUX8(jj,ii) 
	  end do
	end do
	
	return
	END SUBROUTINE  
     	

            
!     ......................................................
      SUBROUTINE F2_0(RO,AKAPA10,AKAPA01,AFME,AF2ME,AF2_MV,NDIM)
!     ......................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
	DIMENSION RO(NDIM,NDIM)
	DIMENSION AKAPA10(NDIM,NDIM)
	DIMENSION AKAPA01(NDIM,NDIM)
	DIMENSION AFME(NDIM,NDIM)
	DIMENSION AF2ME(NDIM,NDIM)

	DIMENSION AUX1(NDIM,NDIM)
	DIMENSION AUX2(NDIM,NDIM)
	DIMENSION AUX3(NDIM,NDIM)
	DIMENSION AUX4(NDIM,NDIM)
	DIMENSION AUX5(NDIM,NDIM)
	DIMENSION AUX6(NDIM,NDIM)


	call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	AF2ME,NDIM,ro,NDIM,zero,AUX1,NDIM)  

	call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	AFME,NDIM,ro,NDIM,zero,AUX2,NDIM)  

	call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	AFME,NDIM,AKAPA10,NDIM,zero,AUX3,NDIM)  

	call DGEMM ('t','t',NDIM,NDIM,NDIM,one,   &
     &	AFME,NDIM,AKAPA01,NDIM,zero,AUX4,NDIM)  

	call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	AUX2,NDIM,AUX2,NDIM,zero,AUX5,NDIM)  

	call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	AUX3,NDIM,AUX4,NDIM,zero,AUX6,NDIM)  

	
	tr1=zero
	tr2=zero
	tr3=zero
	tr4=zero
	do ii=1,NDIM
         tr1=tr1+AUX1(ii,ii)
         tr2=tr2+AUX2(ii,ii)
         tr3=tr3+AUX5(ii,ii)
         tr4=tr4+AUX6(ii,ii)
	end do
	
	AF2_MV=tr1+tr2*tr2-tr3+tr4
	

!       print *, 'tr1=',tr1
!       print *, 'tr2*tr2=',tr2*tr2
!       print *, '-tr3=',-tr3
!       print *, 'tr4=',tr4
     
        end subroutine




