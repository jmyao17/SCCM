        subroutine HFB_ENER_comp(zro,zakapa01,zakin,zgamma,zdelta10,zEkin,zEHF,zEPa,zEHFB,NDIM)
	implicit real*8 (a-h,o-y)
	implicit complex*16 (z)
	
	PARAMETER (zone=(1.d0,0.d0))
	PARAMETER (zzero=(0.d0,0.d0))
	
	DIMENSION zro(NDIM,NDIM)
	DIMENSION zakapa01(NDIM,NDIM)
	DIMENSION zakin(NDIM,NDIM)
	DIMENSION zgamma(NDIM,NDIM)
	DIMENSION zdelta10(NDIM,NDIM)
	DIMENSION zAUX1(NDIM,NDIM)
	DIMENSION zAUX2(NDIM,NDIM)
	DIMENSION zAUX3(NDIM,NDIM)
	
	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,zakin,NDIM,zro,NDIM,zzero,zAUX1,NDIM)  
	call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,zgamma,NDIM,zro,NDIM,zzero,zAUX2,NDIM)  
	call ZGEMM ('n','t',NDIM,NDIM,NDIM,zone,zdelta10,NDIM,zakapa01,NDIM,zzero,zAUX3,NDIM) 
        
	zEHFB=zzero
	zEHF=zzero
	zEPa=zzero
	zEkin=zzero
        do ii=1,NDIM
         zEkin=zEkin+zAUX1(ii,ii)
	 zEHF=zEHF+.5*(zAUX2(ii,ii))
	 zEPa=zEPa+.5*(zAUX3(ii,ii))
	end do 
	zEHFB=zEkin+zEHF+zEPa
    
        end subroutine
	


