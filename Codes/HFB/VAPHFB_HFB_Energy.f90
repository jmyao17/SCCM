!  ................................................
        subroutine HFB_ENER(ro,akapa01,akin,gamma,delta10,Ekin,EHF,EPa,EHFB,NDIM)
!  ................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
        DIMENSION ro(NDIM,NDIM)
        DIMENSION akapa01(NDIM,NDIM)
        DIMENSION akin(NDIM,NDIM)
        DIMENSION gamma(NDIM,NDIM)
        DIMENSION ham(NDIM,NDIM)
        DIMENSION delta10(NDIM,NDIM)
        DIMENSION AUX1(NDIM,NDIM)
        DIMENSION AUX2(NDIM,NDIM)
        DIMENSION AUX3(NDIM,NDIM)
	
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	akin,NDIM,ro,NDIM,zero,AUX1,NDIM)  
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,   &
     &	gamma,NDIM,ro,NDIM,zero,AUX2,NDIM)  
        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,   &
     &	delta10,NDIM,akapa01,NDIM,zero,AUX3,NDIM) 
        
        EHFB=zero
        EHF=zero
        EPa=zero
        Ekin=zero
        do ii=1,NDIM
         Ekin=Ekin+AUX1(ii,ii)
         EHF=EHF+.5*(AUX2(ii,ii))
         EPa=EPa+.5*(AUX3(ii,ii))
        end do 
        EHFB=Ekin+EHF+EPa
    
     
        end subroutine



        subroutine HFB_ENER_COMPLEX(zro,zkapa01,zakin,zgamma,zdelta10,ZEHFB,NDIM)
        use VAPHFB_Par
        implicit real*8 (a-h,o-y)
        implicit complex*16(z)

        DIMENSION zro(NDIM,NDIM)
        DIMENSION zkapa01(NDIM,NDIM)
        DIMENSION zakin(NDIM,NDIM)
        DIMENSION zgamma(NDIM,NDIM)
        DIMENSION zdelta10(NDIM,NDIM)
        DIMENSION ZAUX1(NDIM,NDIM)
        DIMENSION ZAUX2(NDIM,NDIM)
        DIMENSION ZAUX3(NDIM,NDIM)
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,      &
     &	zakin,NDIM,zro,NDIM,zzero,ZAUX1,NDIM)  

        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,      &
     &	zgamma,NDIM,zro,NDIM,zzero,ZAUX2,NDIM)  

        call ZGEMM ('n','t',NDIM,NDIM,NDIM,zone,      &
     &	zdelta10,NDIM,zkapa01,NDIM,zzero,ZAUX3,NDIM) 
        
        ZEHFB=zzero
        ZEHF=zzero
        ZEPa=zzero
        ZEkin=zzero

        do ii=1,NDIM
         ZEkin=ZEkin+ZAUX1(ii,ii)
         ZEHF=ZEHF+.5*(ZAUX2(ii,ii))
         ZEPa=ZEPa+.5*(ZAUX3(ii,ii))
        end do 

        ZEHFB=ZEkin+ZEHF+ZEPa
     
        end subroutine


        subroutine HFB_ENER_fin(ro,akapa01,akin,gamma,delta10,     &
     &  Ekin_P,EHF_P,EPa_P,EHFB_P,                                 &
     &  Ekin_N,EHF_N,EPa_N,EHFB_N,NDIM)

        USE VAPHFB_PAR
        implicit real*8 (a-h,o-z)
	
        DIMENSION ro(NDIM,NDIM)
        DIMENSION ropp(NDIM/2,NDIM/2)
        DIMENSION ropn(NDIM/2,NDIM/2)
        DIMENSION ronp(NDIM/2,NDIM/2)
        DIMENSION ronn(NDIM/2,NDIM/2)
        DIMENSION akapa01(NDIM,NDIM)
        DIMENSION akapa01pp(NDIM/2,NDIM/2)
        DIMENSION akapa01pn(NDIM/2,NDIM/2)
        DIMENSION akapa01np(NDIM/2,NDIM/2)
        DIMENSION akapa01nn(NDIM/2,NDIM/2)
        DIMENSION akin(NDIM,NDIM)
        DIMENSION gamma(NDIM,NDIM)
        DIMENSION gammapp(NDIM/2,NDIM/2)
        DIMENSION gammapn(NDIM/2,NDIM/2)
        DIMENSION gammanp(NDIM/2,NDIM/2)
        DIMENSION gammann(NDIM/2,NDIM/2)
        DIMENSION delta10pp(NDIM/2,NDIM/2)
        DIMENSION delta10pn(NDIM/2,NDIM/2)
        DIMENSION delta10np(NDIM/2,NDIM/2)
        DIMENSION delta10nn(NDIM/2,NDIM/2)
        DIMENSION ham(NDIM,NDIM)
        DIMENSION delta10(NDIM,NDIM)
        DIMENSION AUX1(NDIM,NDIM)
        DIMENSION AUX2(NDIM,NDIM)
        DIMENSION AUX3(NDIM,NDIM)


        DIMENSION AUXHFpp(NDIM/2,NDIM/2)
        DIMENSION AUXHFpn(NDIM/2,NDIM/2)
        DIMENSION AUXHFnp(NDIM/2,NDIM/2)
        DIMENSION AUXHFnn(NDIM/2,NDIM/2)
        DIMENSION AUXPapp(NDIM/2,NDIM/2)
        DIMENSION AUXPapn(NDIM/2,NDIM/2)
        DIMENSION AUXPanp(NDIM/2,NDIM/2)
        DIMENSION AUXPann(NDIM/2,NDIM/2)
        
        ip=0
        in=NDIM/2
        
        do ii=1,NDIM/2
         do jj=1,NDIM/2
           gammapp(ii,jj)=gamma(ii+ip,jj+ip)
           gammapn(ii,jj)=gamma(ii+ip,jj+in)
           gammanp(ii,jj)=gamma(ii+in,jj+ip)
           gammann(ii,jj)=gamma(ii+in,jj+in)

           ropp(ii,jj)=ro(ii+ip,jj+ip)
           ropn(ii,jj)=ro(ii+ip,jj+in)
           ronp(ii,jj)=ro(ii+in,jj+ip)
           ronn(ii,jj)=ro(ii+in,jj+in)

           delta10pp(ii,jj)=delta10(ii+ip,jj+ip)
           delta10pn(ii,jj)=delta10(ii+ip,jj+in)
           delta10np(ii,jj)=delta10(ii+in,jj+ip)
           delta10nn(ii,jj)=delta10(ii+in,jj+in)

           akapa01pp(ii,jj)=akapa01(ii+ip,jj+ip)
           akapa01pn(ii,jj)=akapa01(ii+ip,jj+in)
           akapa01np(ii,jj)=akapa01(ii+in,jj+ip)
           akapa01nn(ii,jj)=akapa01(ii+in,jj+in)


           
         end do 
        end do
        
        call DGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	gammapp,NDIM/2,ropp,NDIM/2,zero,AUXHFpp,NDIM/2)  
        call DGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	gammapn,NDIM/2,ronp,NDIM/2,zero,AUXHFpn,NDIM/2)  

        call DGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	gammanp,NDIM/2,ropn,NDIM/2,zero,AUXHFnp,NDIM/2)  
        call DGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	gammann,NDIM/2,ronn,NDIM/2,zero,AUXHFnn,NDIM/2)  

        call DGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	delta10pp,NDIM/2,akapa01pp,NDIM/2,zero,AUXPapp,NDIM/2)  
        call DGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	delta10pn,NDIM/2,akapa01pn,NDIM/2,zero,AUXPapn,NDIM/2)  

        call DGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	delta10np,NDIM/2,akapa01np,NDIM/2,zero,AUXPanp,NDIM/2)  
        call DGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,one,        &
     &	delta10nn,NDIM/2,akapa01nn,NDIM/2,zero,AUXPann,NDIM/2)  

     
        EHFpp=0.d0
        EHFpn=0.d0
        EHFnp=0.d0
        EHFnn=0.d0
        EPapp=0.d0
        EPapn=0.d0
        EPanp=0.d0
        EPann=0.d0

        do ii=1,NDIM/2
         EHFpp=EHFpp+0.5*AUXHFpp(ii,ii)
         EHFpn=EHFpn+0.5*AUXHFpn(ii,ii)
         EHFnp=EHFnp+0.5*AUXHFnp(ii,ii)
         EHFnn=EHFnn+0.5*AUXHFnn(ii,ii)
         EPapp=EPapp+0.5*AUXPapp(ii,ii)
         EPapn=EPapn+0.5*AUXPapn(ii,ii)
         EPanp=EPanp+0.5*AUXPanp(ii,ii)
         EPann=EPann+0.5*AUXPann(ii,ii)
        end do
	
	
	
	
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,  &
     &	akin,NDIM,ro,NDIM,zero,AUX1,NDIM)  
        call DGEMM ('n','n',NDIM,NDIM,NDIM,one,  &
     &	gamma,NDIM,ro,NDIM,zero,AUX2,NDIM)  
        call DGEMM ('n','t',NDIM,NDIM,NDIM,one,  &
     &	delta10,NDIM,akapa01,NDIM,zero,AUX3,NDIM) 
        
        EHFB_P=zero
        EHF_P=zero
        EPa_P=zero
        Ekin_P=zero
        EHFB_N=zero
        EHF_N=zero
        EPa_N=zero
        Ekin_N=zero
            do ii=1,NDIM/2
             Ekin_P=Ekin_P+AUX1(ii,ii)
         EHF_P=EHF_P+.5*(AUX2(ii,ii))
         EPa_P=EPa_P+.5*(AUX3(ii,ii))
             Ekin_N=Ekin_N+AUX1(ii+NDIM/2,ii+NDIM/2)
         EHF_N=EHF_N+.5*(AUX2(ii+NDIM/2,ii+NDIM/2))
         EPa_N=EPa_N+.5*(AUX3(ii+NDIM/2,ii+NDIM/2))
        end do 
            EHFB_P=Ekin_P+EHF_P+EPa_P
            EHFB_N=Ekin_N+EHF_N+EPa_N
        
        
	write(6,'(A20,5A15)') ' ','prot-prot','neut-neut',  &
     &	'neut-prot','prot-neut','total'
        write(6,'(A20,5A15)') '','------','------',         &
     &	'------','------','------'
        write(6,*) ' '
        write(6,'(A20,5F15.6)') 'Kinetic Energy',Ekin_P,Ekin_N,  &
     &	zero,zero,Ekin_P+Ekin_N

        write(6,'(A20,5F15.6)') 'Hartree-Fock Energy',EHFpp,EHFnn, &
     &	EHFnp,EHFpn,EHFpp+EHFpn+EHFnp+EHFnn	

        write(6,'(A20,5F15.6)') 'Pairing Energy',EPapp,EPann,      &
     &	EPanp,EPapn,EPapp+EPapn+EPanp+EPann
        write(6,*) ' '
        end subroutine



        subroutine HFB_ENER_COMPLEX_fin(zro,zkapa01,zakin,  &
     &  zgamma,zdelta10,ZEkin_P,ZEkin_N,                    &
     &  ZEHFPP,ZEHFPN,ZEHFNP,ZEHFNN,                        &
     &  ZEPaPP,ZEPaPN,ZEPaNP,ZEPaNN,                        &
     &  ZEHFB_P,ZEHFB_N,NDIM)

        Use VAPHFB_Par
        implicit real*8 (a-h,o-y)
	implicit complex*16 (z)
	
	
        DIMENSION zro(NDIM,NDIM)
        DIMENSION zropp(NDIM/2,NDIM/2)
        DIMENSION zropn(NDIM/2,NDIM/2)
        DIMENSION zronp(NDIM/2,NDIM/2)
        DIMENSION zronn(NDIM/2,NDIM/2)
        DIMENSION zkapa01(NDIM,NDIM)
        DIMENSION zkapa01pp(NDIM/2,NDIM/2)
        DIMENSION zkapa01pn(NDIM/2,NDIM/2)
        DIMENSION zkapa01np(NDIM/2,NDIM/2)
        DIMENSION zkapa01nn(NDIM/2,NDIM/2)
        DIMENSION zakin(NDIM,NDIM)
        DIMENSION zgamma(NDIM,NDIM)
        DIMENSION zgammapp(NDIM/2,NDIM/2)
        DIMENSION zgammapn(NDIM/2,NDIM/2)
        DIMENSION zgammanp(NDIM/2,NDIM/2)
        DIMENSION zgammann(NDIM/2,NDIM/2)
        DIMENSION zdelta10pp(NDIM/2,NDIM/2)
        DIMENSION zdelta10pn(NDIM/2,NDIM/2)
        DIMENSION zdelta10np(NDIM/2,NDIM/2)
        DIMENSION zdelta10nn(NDIM/2,NDIM/2)
        DIMENSION zham(NDIM,NDIM)
        DIMENSION zdelta10(NDIM,NDIM)
        DIMENSION zAUX1(NDIM,NDIM)
        DIMENSION zAUX2(NDIM,NDIM)
        DIMENSION zAUX3(NDIM,NDIM)


        DIMENSION zAUXHFpp(NDIM/2,NDIM/2)
        DIMENSION zAUXHFpn(NDIM/2,NDIM/2)
        DIMENSION zAUXHFnp(NDIM/2,NDIM/2)
        DIMENSION zAUXHFnn(NDIM/2,NDIM/2)
        DIMENSION zAUXPapp(NDIM/2,NDIM/2)
        DIMENSION zAUXPapn(NDIM/2,NDIM/2)
        DIMENSION zAUXPanp(NDIM/2,NDIM/2)
        DIMENSION zAUXPann(NDIM/2,NDIM/2)
        
        ip=0
        in=NDIM/2
        
        do ii=1,NDIM/2
         do jj=1,NDIM/2
           zgammapp(ii,jj)=zgamma(ii+ip,jj+ip)
           zgammapn(ii,jj)=zgamma(ii+ip,jj+in)
           zgammanp(ii,jj)=zgamma(ii+in,jj+ip)
           zgammann(ii,jj)=zgamma(ii+in,jj+in)

           zropp(ii,jj)=zro(ii+ip,jj+ip)
           zropn(ii,jj)=zro(ii+ip,jj+in)
           zronp(ii,jj)=zro(ii+in,jj+ip)
           zronn(ii,jj)=zro(ii+in,jj+in)

           zdelta10pp(ii,jj)=zdelta10(ii+ip,jj+ip)
           zdelta10pn(ii,jj)=zdelta10(ii+ip,jj+in)
           zdelta10np(ii,jj)=zdelta10(ii+in,jj+ip)
           zdelta10nn(ii,jj)=zdelta10(ii+in,jj+in)

           zkapa01pp(ii,jj)=zkapa01(ii+ip,jj+ip)
           zkapa01pn(ii,jj)=zkapa01(ii+ip,jj+in)
           zkapa01np(ii,jj)=zkapa01(ii+in,jj+ip)
           zkapa01nn(ii,jj)=zkapa01(ii+in,jj+in)


           
         end do 
        end do
        
        call ZGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zgammapp,NDIM/2,zropp,NDIM/2,zzero,ZAUXHFpp,NDIM/2)  
        call ZGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zgammapn,NDIM/2,zronp,NDIM/2,zzero,ZAUXHFpn,NDIM/2)  

        call ZGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zgammanp,NDIM/2,zropn,NDIM/2,zzero,ZAUXHFnp,NDIM/2)  
        call ZGEMM ('n','n',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zgammann,NDIM/2,zronn,NDIM/2,zzero,ZAUXHFnn,NDIM/2)  

        call ZGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zdelta10pp,NDIM/2,zkapa01pp,NDIM/2,zzero,ZAUXPapp,NDIM/2)  
        call ZGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zdelta10pn,NDIM/2,zkapa01pn,NDIM/2,zzero,ZAUXPapn,NDIM/2)  

        call ZGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zdelta10np,NDIM/2,zkapa01np,NDIM/2,zzero,ZAUXPanp,NDIM/2)  
        call ZGEMM ('n','t',NDIM/2,NDIM/2,NDIM/2,zone,        &
     &	zdelta10nn,NDIM/2,zkapa01nn,NDIM/2,zzero,ZAUXPann,NDIM/2)  

     
        ZEHFpp=zzero
        ZEHFpn=zzero
        ZEHFnp=zzero
        ZEHFnn=zzero
        ZEPapp=zzero
        ZEPapn=zzero
        ZEPanp=zzero
        ZEPann=zzero

        do ii=1,NDIM/2
         ZEHFpp=ZEHFpp+0.5*ZAUXHFpp(ii,ii)
         ZEHFpn=ZEHFpn+0.5*ZAUXHFpn(ii,ii)
         ZEHFnp=ZEHFnp+0.5*ZAUXHFnp(ii,ii)
         ZEHFnn=ZEHFnn+0.5*ZAUXHFnn(ii,ii)
         ZEPapp=ZEPapp+0.5*ZAUXPapp(ii,ii)
         ZEPapn=ZEPapn+0.5*ZAUXPapn(ii,ii)
         ZEPanp=ZEPanp+0.5*ZAUXPanp(ii,ii)
         ZEPann=ZEPann+0.5*ZAUXPann(ii,ii)
        end do
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,         &
     &	zakin,NDIM,zro,NDIM,zzero,ZAUX1,NDIM)  
        call ZGEMM ('n','n',NDIM,NDIM,NDIM,zone,         &
     &	zgamma,NDIM,zro,NDIM,zzero,ZAUX2,NDIM)  
        call ZGEMM ('n','t',NDIM,NDIM,NDIM,zone,         &
     &	zdelta10,NDIM,zkapa01,NDIM,zzero,ZAUX3,NDIM) 
        
        ZEHFB_P=zzero
        ZEHF_P=zzero
        ZEPa_P=zzero
        ZEkin_P=zzero
        ZEHFB_N=zzero
        ZEHF_N=zzero
        ZEPa_N=zzero
        ZEkin_N=zzero
            do ii=1,NDIM/2
             ZEkin_P=ZEkin_P+ZAUX1(ii,ii)
         ZEHF_P=ZEHF_P+.5*(ZAUX2(ii,ii))
         ZEPa_P=ZEPa_P+.5*(ZAUX3(ii,ii))
             ZEkin_N=ZEkin_N+ZAUX1(ii+NDIM/2,ii+NDIM/2)
         ZEHF_N=ZEHF_N+.5*(ZAUX2(ii+NDIM/2,ii+NDIM/2))
         ZEPa_N=ZEPa_N+.5*(ZAUX3(ii+NDIM/2,ii+NDIM/2))
        end do 
            ZEHFB_P=ZEkin_P+ZEHF_P+ZEPa_P
            ZEHFB_N=ZEkin_N+ZEHF_N+ZEPa_N
        
        end subroutine

