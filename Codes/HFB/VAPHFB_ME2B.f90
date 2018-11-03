        SUBROUTINE G_20(U0,V0,AGME,AG_20,NLEV)
        USE VAPHFB_PAR 
 
        implicit real*8 (a-h,o-z)
        
        DIMENSION V0(NLEV,NLEV)
        DIMENSION U0(NLEV,NLEV)
        DIMENSION AGME(NLEV,NLEV)
        DIMENSION AG_20(NLEV,NLEV)
        DIMENSION AG_20AUX(NLEV,NLEV)
        DIMENSION AUX1(NLEV,NLEV),AUX2(NLEV,NLEV)
        DIMENSION AUX3(NLEV,NLEV),AUX4(NLEV,NLEV)
        
        
        call DGEMM ('t','n',NLEV,NLEV,NLEV,one,U0,NLEV,AGME,NLEV,zero,AUX1,NLEV)       
        call DGEMM ('n','n',NLEV,NLEV,NLEV,one,AUX1,NLEV,U0,NLEV,zero,AUX3,NLEV)       

        call DGEMM ('t','n',NLEV,NLEV,NLEV,one,V0,NLEV,AGME,NLEV,zero,AUX2,NLEV)       
        call DGEMM ('n','n',NLEV,NLEV,NLEV,one,AUX2,NLEV,V0,NLEV,zero,AUX4,NLEV)

        do ii=1,NLEV
         do jj=1,NLEV
          AG_20AUX(jj,ii)=AUX3(jj,ii)-AUX4(jj,ii)
         end do 
        end do
        
        do ii=1,NLEV
         do jj=1,NLEV
          AG_20(jj,ii)=0.5*(AG_20AUX(jj,ii)  -AG_20AUX(ii,jj))     
         end do 
        end do


        END SUBROUTINE
! ....................................... 
       SUBROUTINE F_20(U0,V0,AFME,AF_20,NLEV)
! ....................................... 
        USE VAPHFB_PAR 
        implicit real*8 (a-h,o-z)

        DIMENSION V0(NLEV,NLEV)
        DIMENSION U0(NLEV,NLEV)
        DIMENSION AFME(NLEV,NLEV)
        DIMENSION AF_20(NLEV,NLEV)
        DIMENSION AF_20AUX(NLEV,NLEV)
        DIMENSION AUX1(NLEV,NLEV),AUX2(NLEV,NLEV)
        DIMENSION AUX3(NLEV,NLEV),AUX4(NLEV,NLEV)


        call DGEMM ('t','n',NLEV,NLEV,NLEV,one,          &
     &  U0,NLEV,AFME,NLEV,zero,AUX1,NLEV)
        call DGEMM ('n','n',NLEV,NLEV,NLEV,one,          &
     &  AUX1,NLEV,V0,NLEV,zero,AUX3,NLEV)     

        call DGEMM ('t','t',NLEV,NLEV,NLEV,one,          &
     &  V0,NLEV,AFME,NLEV,zero,AUX2,NLEV)      
        call DGEMM ('n','n',NLEV,NLEV,NLEV,one,          &
     &  AUX2,NLEV,U0,NLEV,zero,AUX4,NLEV)
     
        do ii=1,NLEV
	 do jj=1,NLEV
	  AF_20AUX(jj,ii)=AUX3(jj,ii)-AUX4(jj,ii)
	 end do
	end do

        do ii=1,NLEV
	 do jj=1,NLEV
	  AF_20(jj,ii)=0.5*(AF_20AUX(jj,ii) -AF_20AUX(ii,jj)) 
	 end do
	end do
     
        END SUBROUTINE 
