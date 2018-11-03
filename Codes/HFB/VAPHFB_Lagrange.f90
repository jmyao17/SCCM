
!      .................................................................
       SUBROUTINE LAGRANGE_MULT(NCONS,H20,ACON20,alagCON_0,alagCON,NDIM)
!      .................................................................
!      INPUTS: NCONS,H20,ACON20
!      .................................................................
       USE VAPHFB_PAR 
       implicit real*8 (a-h,o-z)
        
        DIMENSION alagCON_0(NCONSMAX)
        DIMENSION alagCON(NCONSMAX)
!       ...................................
        DIMENSION ACON20(NCONS*NDIM*NDIM)
        DIMENSION H20(NDIM,NDIM)
        DIMENSION AUX1(((NCONS+1)*NCONS/2)*NDIM*NDIM)
        DIMENSION AUX2(NCONS*NDIM*NDIM)
        DIMENSION Amatrix(NCONS,NCONS),Bmatrix(NCONS)
        DIMENSION IPIV(NCONS)
           
        ll=0
        do ii=1,NCONS
         indi=(ii-1)*NDIM*NDIM+1

!        ............ number of constraint terms
         do jj=ii,NCONS
            indj=(jj-1)*NDIM*NDIM+1
            ll=ll+1
            indl=(ll-1)*NDIM*NDIM+1

!      ............ AUX1 = ACON20 * ACON20
        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,ACON20(indi), &
     &               NDIM,ACON20(indj),NDIM,0.d0,AUX1(indl),NDIM) 
          
         end do

!      ............ AUX2 = H20 * ACON20
         call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,               &
     &        H20,NDIM,ACON20(indi),NDIM,0.d0,AUX2(indi),NDIM)

        end do
        
        
        ll=0
        do ii=1,NCONS
         indi=(ii-1)*NDIM*NDIM+1
         do jj=ii,NCONS
         ll=ll+1
         indl=(ll-1)*NDIM*NDIM+1

           Amatrix(jj,ii)=zero   
           do mm=1,NDIM
              idig=(mm-1)*NDIM+mm-1
              Amatrix(jj,ii)=Amatrix(jj,ii)+AUX1(indl+idig)   
           end do
              Amatrix(ii,jj)=Amatrix(jj,ii)           
         end do

         Bmatrix(ii)=zero
           do mm=1,NDIM
              idig=(mm-1)*NDIM+mm-1
              Bmatrix(ii)=Bmatrix(ii)+AUX2(indi+idig)	   
           end do
         
        end do


        call DGESV(NCONS,1,Amatrix,NCONS,IPIV,Bmatrix,NCONS,INFO)
            
        do ii=1,NCONS
             alagCON(ii)=Bmatrix(ii)
        end do
        
        if (INFO.ne.0) then
         do ii=1,NCONS
          alagCON(ii)=alagCON_0(ii)
         end do
        end if

        do ii=1,NCONS
          alagCON_0(ii)=alagCON(ii)
        end do

        
        END SUBROUTINE  








         
