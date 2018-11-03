!       .....................................
        SUBROUTINE N_0(V0,AN0,AZ0,NDIM)  
!       .....................................
!       Particle Number
!       .....................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        DIMENSION V0(NDIM,NDIM)
        DIMENSION ro(NDIM,NDIM)
        
        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,     &
     &              V0,NDIM,V0,NDIM,0.d0,ro,NDIM)
        AZ0=zero
        AN0=zero
        do ii=1,NDIM/2
           AZ0=AZ0+ro(ii,ii)
           AN0=AN0+ro(ii+NDIM/2,ii+NDIM/2)
        end do
       END SUBROUTINE 

!       .....................................
        SUBROUTINE betgamm_0(V0,bet2t,gam2t,NDIM)
!       .....................................
!       deformation: beta2,gamma2 
!       .....................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        DIMENSION V0(NDIM,NDIM)
        DIMENSION ro(NDIM,NDIM)

        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,     &
     &              V0,NDIM,V0,NDIM,0.d0,ro,NDIM)
!      ........ deformation parameter
        call Q_0(cME1B%Q20t,ro,Q20p,Q20n,NDIM)
        call Q_0(cME1B%Q2_2t,ro,Q2_2p,Q2_2n,NDIM)
        call Q_0(cME1B%Q22t,ro,Q22p,Q22n,NDIM)
        Q20tbeta  = Q20p+Q20n
        Q2_2tbeta = Q2_2p+Q2_2n
        Q22tbeta  = Q22p+Q22n
!     ............. determination of beta2 from Q20 and Q22
        bsignt = Q20tbeta/abs(Q20tbeta)
        bet2t=bsignt*Const%Q2BA * dsqrt(Q20tbeta**2+Q22tbeta**2+Q2_2tbeta**2)
        gam2t  = 180.d0*atan(sqrt(2.d0)*Q22tbeta/Q20tbeta)/pi
        return
        end
!   ..............................................................
        SUBROUTINE UV2Density(U0,V0,RO,akapa10,akapa01,NDIM)  
!   ..............................................................
!    Densities rho and kappa
!    rho = V*V^T; kappa_10= V*U^T; kappa_01 = U*V^T  
!   ..............................................................
	USE VAPHFB_PAR
	implicit real*8 (a-h,o-z)
	
        DIMENSION V0(NDIM,NDIM)
        DIMENSION U0(NDIM,NDIM)
	DIMENSION ro(NDIM,NDIM)
	DIMENSION akapa10(NDIM,NDIM),akapa01(NDIM,NDIM)

        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,       &
     &             V0,NDIM,V0,NDIM,0.d0,ro,NDIM)
 
        call DGEMM ('n','t',NDIM,NDIM,NDIM,1.d0,       &
     &  V0,NDIM,U0,NDIM,0.d0,akapa10,NDIM)  
       
        call DGEMM ('n','t',NDIM,NDIM,NDIM,-1.d0,      &
     &  U0,NDIM,V0,NDIM,0.d0,akapa01,NDIM)

      if(Input%NPMix .eq. 0) then
      do ii=1,NDIM
      do jj=1,NDIM
!  ................w/o np-mixing and pairing
         if(tnljm%t(ii) .ne. tnljm%t(jj)) then
            ro(ii,jj)       = 0.d0
            akapa10(ii,jj)  = 0.d0
            akapa01(ii,jj)  = 0.d0
          endif
       enddo
       enddo
       endif

        END SUBROUTINE 
            

!   ..............................................................
      SUBROUTINE Q_0(QME,RO,Q_0p,Q_0n,NDIM)  
!   ..............................................................
!   calculate the expectation value of 1B operator:
!   Q_0p = sum_ij QME(i,j)*Rho_p(j,i)
!   Q_0n = sum_ij QME(i,j)*Rho_n(j,i)
!   ..............................................................
        USE VAPHFB_PAR
	implicit real*8 (a-h,o-z)
	
        DIMENSION QME(NDIM,NDIM)
	DIMENSION RO(NDIM,NDIM)
	DIMENSION AUX(NDIM,NDIM)
        
	
	call DGEMM ('n','n',NDIM,NDIM,NDIM,1.d0,        &
     &	QME,NDIM,RO,NDIM,0.d0,AUX,NDIM)	
	
       	
	Q_0p=zero
	Q_0n=zero
	do ii=1,NDIM/2
	 Q_0p=Q_0p+AUX(ii,ii)
	 Q_0n=Q_0n+AUX(ii+NDIM/2,ii+NDIM/2)
	end do
	
       	END SUBROUTINE 


      SUBROUTINE ZN_0(ZROBP,Zprot_0,Zneut_0,N)
      USE VAPHFB_PAR
      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
!      Parameter (Zone=(1.d0,0.d0))
!      Parameter (Zzero=(0.d0,0.d0))
!      Parameter (Zimag=(0.d0,1.d0))
      
      Dimension ZROBP(N,N)
      
      zprot_0=zzero
      zneut_0=zzero
      
      do ii=1,N/2
       zprot_0=zprot_0+ZROBP(ii,ii)
       zneut_0=zneut_0+ZROBP(ii+N/2,ii+N/2) 
      end do
      
      
      end subroutine


!      .................................................
        SUBROUTINE ZQ_0(ZQME,ZRO,zQ_0p,zQ_0n,NDIM)  
!      .................................................
        USE VAPHFB_PAR        
        implicit real*8 (a-h,o-y)
        implicit complex*16 (z)
!        Parameter (zzero=(0.d0,0.d0))
!        Parameter (zone=(1.d0,0.d0))
        
        DIMENSION ZQME(NDIM,NDIM)
        DIMENSION ZRO(NDIM,NDIM)
        DIMENSION ZAUX(NDIM,NDIM)
            
        
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,     &
     &       ZQME,NDIM,ZRO,NDIM,zzero,ZAUX,NDIM)
        zQ_0p=zzero
        zQ_0n=zzero
        do ii=1,NDIM/2
         zQ_0p=zQ_0p+ZAUX(ii,ii)
         zQ_0n=zQ_0n+ZAUX(ii+NDIM/2,ii+NDIM/2)
        end do
        END SUBROUTINE 
