! ......................................................................
        subroutine NZPROJ_VAP(AZ_P,AN_P,EPNP,H20_P,NFOM,NLEV)
! ......................................................................
        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16(z)
        
        DIMENSION H20_P(NLEV,NLEV)

        integer LLPN

        DIMENSION ZU0(NLEV,NLEV)
        DIMENSION ZV0(NLEV,NLEV)
        complex*16, dimension (:,:), allocatable :: ZU0bar,ZV0bar,ZU0tilde,ZV0tilde
        complex*16, dimension (:,:), allocatable :: zgamma,zdelta10,zdelta01,zham
        DIMENSION ZZc(NLEV,NLEV)
        DIMENSION ZRO(NLEV,NLEV)
        DIMENSION Zkapa10(NLEV,NLEV)
        DIMENSION Zkapa01(NLEV,NLEV)
!       ...gradients
        DIMENSION ZH20phi(NLEV,NLEV)
        DIMENSION Zsum_Aphi(NLEV,NLEV)
        DIMENSION Zsum_H20phi(NLEV,NLEV)
        DIMENSION Zsum_Aphi_H20phi(NLEV,NLEV)
        DIMENSION Zsum_AphiA(NLEV,NLEV)
        DIMENSION Zsum_H20phiA(NLEV,NLEV)
        DIMENSION Zsum_Aphi_H20phiA(NLEV,NLEV)
        DIMENSION ZAphi(NLEV,NLEV)
        DIMENSION Zakin(NLEV,NLEV)
        dimension Nphip(1:NFOM*NFOM)
        dimension Nphin(1:NFOM*NFOM)
        Dimension Zsumphi(NOBS)
        

         Aprot = Input%nprot 
         Aneut = Input%nneut 
        
        do ii=1,NOBS
           Zsumphi(ii)=zzero
        end do
        
        do ii=1,NLEV
         do jj=1,NLEV
            ZV0(jj,ii)=zone*HFB%V0(jj,ii)
            ZU0(jj,ii)=zone*HFB%U0(jj,ii)
            ZAkin(jj,ii)=zone*H%ME1BM(jj,ii)
         end do
        end do

        zDetvac=zone
        if(NFOM .gt. 1) then
           call zbareover(ZV0,ZU0,ZZc,Det,NLEV)
           zDetvac=dsqrt((1./Det))*zone
         endif
        aaa=1.d0*NLEV*(NLEV+1)/2.d0
        snpfaf=(-1.d0)**(aaa)

!!!!!!! initializing gradients
        do ii=1,NLEV
          do jj=1,NLEV
	     Zsum_Aphi(jj,ii)       =zzero
	     Zsum_H20phi(jj,ii)     =zzero
	     Zsum_Aphi_H20phi(jj,ii)=zzero
	     H20_P(jj,ii)           =0.d0
	 end do
	end do
        
        
      if(.NOT. ALLOCATED(ZU0tilde))  ALLOCATE(ZU0tilde(NLEV,NLEV))
      if(.NOT. ALLOCATED(ZV0tilde))  ALLOCATE(ZV0tilde(NLEV,NLEV))
      if(.NOT. ALLOCATED(ZU0bar))  ALLOCATE(ZU0bar(NLEV,NLEV))
      if(.NOT. ALLOCATED(ZV0bar))  ALLOCATE(ZV0bar(NLEV,NLEV))
                     
      if(.NOT. ALLOCATED(Zgamma))    ALLOCATE(Zgamma(NLEV,NLEV))
      if(.NOT. ALLOCATED(Zdelta01))  ALLOCATE(Zdelta01(NLEV,NLEV))
      if(.NOT. ALLOCATED(Zdelta10))  ALLOCATE(Zdelta10(NLEV,NLEV))
      if(.NOT. ALLOCATED(Zham))      ALLOCATE(Zham(NLEV,NLEV))


           LLPN=0
        do LLLP=0,NFOM-1
        do LLLN=0,NFOM-1
           LLPN=LLPN+1
           Nphip(LLPN) = LLLP
           Nphin(LLPN) = LLLN
        enddo ! LLPN
        enddo ! LLPN
       DO LLPN=1,NFOM**2   !Fomenko                

        phiLP=(2*pi/NFOM)*Nphip(LLPN)
        if(NFOM.eq.1) phiLP=0.d0*pi/180.d0

        zphas_p=cdexp(-Zimag*phiLP*Aprot)
        phiLN=(2*pi/NFOM)*Nphin(LLPN)
        if(NFOM.eq.1) phiLN=0.d0*pi/180.d0

        zphas_n=cdexp(-Zimag*phiLN*Aneut)

        call UVphi(phiLP,phiLN,ZU0,ZV0,ZU0bar,ZV0bar,ZAphi,ZU0tilde,        &
     &             ZV0tilde,NLEV)



!    norm overlap
      if(NFOM .gt. 1) then
      ! write(*,*) '.... NZPROJ_VAP: Pfaffian ........'      
         call Pfaf_Over(ZZc,ZU0bar,ZV0bar,ZdetPfaf,NLEV)
         znorm =snpfaf*zDetvac*ZdetPfaf*zphas_p*zphas_n
      elseif(NFOM .eq. 1) then
         znorm = zone
      else
        stop 'NFOM should be a positive integer !!'
      endif
        Zsumphi(1)=Zsumphi(1)+znorm

        
        !ro,kappa...
        call ROKAPPA_NZPROJ(ZU0,ZV0,ZU0tilde,ZV0tilde,ZRO,ZKAPA10,ZKAPA01,NLEV)
        
        
        !fields...
!       write(*,*) 'HFB_Field'
       call HFB_FIELD_COMPLEX(zro,zkapa10,ZKAPA01,zakin,zgamma,zham,zdelta10,zdelta01,NLEV)  
!       call HFB_FIELD_COMPLEX_Sph(zro,zkapa10,ZKAPA01,zakin,zgamma,zham,zdelta10,zdelta01,NLEV)  
     
!    ..................... store the gamma   
        if(NFOM.eq.1) then
           HFB%gamma_0 = dreal(zgamma)
           HFB%delta10 = dreal(zdelta10)
           HFB%ham_0   = dreal(zham)
        endif 
       !energy
!        write(*,*) 'HFB_Energy'
       call HFB_ENER_COMPLEX(zro,zkapa01,zakin,zgamma,zdelta10,zEPNP,NLEV)

        zsumphi(2)=zsumphi(2)+zEPNP*znorm
        
        !number of particles for testing purposes
        call ZN_0(ZRO,Zprot_0,Zneut_0,NLEV)

        zsumphi(3)=zsumphi(3)+Zprot_0*znorm
        zsumphi(4)=zsumphi(4)+Zneut_0*znorm
!
!       write(*,*) '  HFB_Grad_phi ...'
       call HFB_GRAD_PHI(ZU0tilde,ZV0tilde,zham,zdelta10,zdelta01,ZH20phi,NLEV)
       
       do ii=1,NLEV
       do jj=1,NLEV
          Zsum_Aphi(jj,ii)  = Zsum_Aphi(jj,ii)+ZAphi(jj,ii)*znorm
          Zsum_H20phi(jj,ii)= Zsum_H20phi(jj,ii)+ZH20phi(jj,ii)*znorm
          Zsum_Aphi_H20phi(jj,ii)= Zsum_Aphi_H20phi(jj,ii)+ZAphi(jj,ii)*zEPNP*znorm

       end do
       end do


        END DO !Fomenko
!    .........................................


!       write(*,*) '  DeAllocate Memory ...'
        DeALLOCATE(ZU0tilde)
        DeALLOCATE(ZV0tilde)
        DeALLOCATE(ZU0bar)
        DeALLOCATE(ZV0bar)
        DeALLOCATE(Zgamma)
        DeALLOCATE(Zdelta10)
        DeALLOCATE(Zdelta01)
        DeALLOCATE(Zham)

       do ii=1,NLEV
       do jj=1,NLEV
	Zsum_AphiA(jj,ii)=0.5*(Zsum_Aphi(jj,ii)-          &
     &                    Zsum_Aphi(ii,jj))/zsumphi(1)

	Zsum_H20phiA(jj,ii)=0.5*(Zsum_H20phi(jj,ii)-      &
     &     Zsum_H20phi(ii,jj))/zsumphi(1)

	Zsum_Aphi_H20phiA(jj,ii)=0.5*(Zsum_Aphi_H20phi(jj,ii)-  &
     &     Zsum_Aphi_H20phi(ii,jj))/zsumphi(1)

       end do
       end do
       EPNP=dreal(zsumphi(2)/zsumphi(1)) 
       AZ_P=dreal(zsumphi(3)/zsumphi(1))   
       AN_P=dreal(zsumphi(4)/zsumphi(1))   

!    .........................................
       tolll=1.e-8
       do ii=1,NLEV
       do jj=1,NLEV
        H20_P(jj,ii)=dreal( Zsum_Aphi_H20phiA(jj,ii)-         &
     &     zsumphi(2)/zsumphi(1)*Zsum_AphiA(jj,ii)+            &
     &     Zsum_H20phiA(jj,ii))
     
        aaaa=dimag(   ( Zsum_Aphi_H20phi(jj,ii) &
     &                 -zsumphi(2)/zsumphi(1)*Zsum_Aphi(jj,ii)  &
     &                 +Zsum_H20phi(jj,ii)                      &
     &                  ) /zsumphi(1) )   
    
         if(aaaa.gt.tolll) then
            write(*,*) 'aaaa= ',aaaa,' tolll=',tolll
!            print*,aaaa  ; stop
         end if  
       end do
       end do
!    .........................................

        end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    !!!!!!!!! Norm ovelap with the bare vacuum !!!!!!!		
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
        

          Subroutine zbareover(ZV,ZU,ZZc,Det,N)
          USE VAPHFB_PAR     
          implicit real*8 (a-h,o-y)
          implicit complex*16 (z)
          
          Dimension ZV(N,N),ZU(N,N)
          Dimension ZUINV(N,N),ZZ(N,N),ZZc(N,N)
          Dimension ZZTZ(N,N)
          Dimension IPIV(N)
          Dimension ZWORK(N)
          Dimension ZWR(N)
          Dimension ZVL(N,N)
          Dimension ZVR(N,N)
          Dimension ZWORK1(2*N)
          Dimension WORK1(2*N)
          
          do ii=1,N
           do jj=1,N
              ZUINV(jj,ii)=ZU(jj,ii)
           end do
          end do

          !LU factorization
          call ZGETRF(N,N,ZUINV,N,IPIV,INFO1)


          if(info1.ne.0) then
          write(6,*) ' In bareover got INFO1 = ',info1,' from ZGETRF'
          stop
          end if
         
          !Inverse
          call ZGETRI(N,ZUINV,N,IPIV,ZWORK,N,INFO2)
          if(info2.ne.0) then
          write(6,*) ' In bareover got INFO2 = ',info2,' from ZGETRI'
          stop
          end if
          

          !Z
          call ZGEMM ('n','n',N,N,N,zone, ZV,N,ZUINV,N,zzero,ZZ,N)
         
          do ii=1,N      
          do jj=1,N
          ZZc(jj,ii)=dconjg(ZZ(jj,ii))
          end do
          end do
          
          
          !ZtZ
          call ZGEMM ('t','n',N,N,N,zone, ZZc,N,ZZc,N,zzero,ZZTZ,N)
                
          !Eigensystem
          call ZGEEV('n','n',N,ZZTZ,N,ZWR,ZVL,N,ZVR,        &
     &                  N,ZWORK1,2*N,WORK1,INFO3 )
          if(info3.ne.0) then
        write(6,*) ' In bareover got INFO3 = ',info3,' from ZGEEV'
        stop
          end if
          
!      do ii=1,N
!       write(20,'(100E20.5)') (Z(ii,jj),jj=1,N)
!       print*,ii,Z(ii,ii)
!      end do

          zdet=zone
          do ii=1,N
           zdet=zdet*(zone+ZWR(ii))
           write(985,*) ii,ZWR(ii)
          end do

          Det=dreal(zdet) !! 1/(<Phi|->^2)
          return
          end  




! .............................................................................
      subroutine UVphi(philZ,philN,ZU,ZV,ZUbar,ZVbar,ZAphi,ZUtilde,ZVtilde,N)
! .............................................................................
      USE VAPHFB_PAR 
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      
      Dimension ZV(N,N),ZU(N,N)
      Dimension ZUbar(N,N),ZVbar(N,N)
      Dimension ZUphi(N,N),ZVphi(N,N)
      Dimension ZUtilde(N,N),ZVtilde(N,N)
      Dimension ZUphi_inv(N,N),ZAphi(N,N),ZAphic(N,N)
      Dimension ZAUX1(N,N),ZAUX2(N,N)

      Dimension IPIV(N)
      Dimension ZWORK(N)

      zphasUZ=cdexp(zimag*philZ)
      zphasUN=cdexp(zimag*philN)
      zphasVZ=cdexp(-zimag*philZ)
      zphasVN=cdexp(-zimag*philN)
      
      do ii=1,N
       do jj=1,N
        ZUbar(jj,ii)=zzero
        ZVbar(jj,ii)=zzero
        end do
      end do
      
      do ii=1,N/2
       do jj=1,N
        ZUbar(ii    ,jj)=zphasUZ*ZU(ii    ,jj)
        ZUbar(ii+N/2,jj)=zphasUN*ZU(ii+N/2,jj)
        ZVbar(ii    ,jj)=zphasVZ*ZV(ii    ,jj)
        ZVbar(ii+N/2,jj)=zphasVN*ZV(ii+N/2,jj)
       end do
      end do

      !Uphi
      call ZGEMM ('t','n',N,N,N,zone,ZV,N,ZVbar,N,zzero,ZUphi,N)

      call ZGEMM ('t','n',N,N,N,zone,ZU,N,ZUbar,N,zone,ZUphi,N)

      !Vphi
      call ZGEMM ('t','n',N,N,N,zone,ZV,N,ZUbar,N,zzero,ZVphi,N)

      call ZGEMM ('t','n',N,N,N,zone,ZU,N,ZVbar,N,zone,ZVphi,N)
     
      !Uphi-1
      do ii=1,N
       do jj=1,N
        ZUphi_inv(jj,ii)=ZUphi(jj,ii)
       end do
      end do
      call ZGETRF(N,N,ZUphi_inv,N,IPIV,INFO1)

      call ZGETRI(N,ZUphi_inv,N,IPIV,ZWORK,N,INFO2)
     
      ! Aphi^*
      call ZGEMM ('n','n',N,N,N,zone,ZVphi,N,ZUphi_inv,N,zzero,ZAphic,N)

      do ii=1,N
       do jj=1,N
        ZAphi(jj,ii)=dconjg(ZAphic(jj,ii))
       end do
      end do


	!Utilde(phi),Vtilde(phi)
      call ZGEMM ('n','n',N,N,N,zone,ZV,N,ZAphi,N,zzero,ZAUX1,N)
      call ZGEMM ('n','n',N,N,N,zone,ZU,N,ZAphi,N,zzero,ZAUX2,N)
     
      do ii=1,N
       do jj=1,N
        ZUtilde(jj,ii)=ZU(jj,ii)+ZAUX1(jj,ii)
	ZVtilde(jj,ii)=ZV(jj,ii)+ZAUX2(jj,ii)
       end do
      end do
      

      
      end subroutine

!     ...............................................................
      subroutine ROKAPPA_NZPROJ(ZU,ZV,ZUtilde,ZVtilde,          &
     & ZRO,ZKAPA10,ZKAPA01,N)
!     ...............................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      
      Dimension ZV(N,N),ZU(N,N)
      Dimension ZUtilde(N,N),ZVtilde(N,N)
      Dimension ZRO(N,N),ZKAPA10(N,N),ZKAPA01(N,N)

     

      call ZGEMM ('n','t',N,N,N,Zone,ZVtilde,N,ZV,N,zzero,ZRO,N)

      call ZGEMM ('n','t',N,N,N,Zone,ZVtilde,N,ZU,N,zzero,ZKAPA10,N)

      call ZGEMM ('n','t',N,N,N,zone,ZV,N,ZUtilde,N,Zzero,ZKAPA01,N) 
     
      
      end subroutine
          







      SUBROUTINE Pfaf_Over(ZM0,ZUbar,ZVbar,ZdetPfaf,N)    
      USE VAPHFB_PAR
      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
      Dimension ZPfaf(2*N,2*N) !rotation matrix 
      Dimension ZUbar(N,N),ZVbar(N,N) 
      Dimension ZUbar_inv(N,N),ZVbar_c(N,N) 
      Dimension ZM0(N,N)
      
      Dimension IPIV(N)
      Dimension ZWORK(N)

      Dimension ZM1(N,N) 
     
      Dimension ZAUX1(N,N) 
      Dimension ZAUX2(N,N) 
      
      DIMENSION Ipiv2(2*N,2)
      
      do ii=1,N
       do jj=1,N
        ZUbar_inv(jj,ii)=dconjg(ZUbar(jj,ii))
	ZVbar_c(jj,ii)=dconjg(ZVbar(jj,ii))
       end do
      end do
      !LU factorization
      call ZGETRF(N,N,ZUbar_inv,N,IPIV,INFO1)


      if(info1.ne.0) then
	write(6,*) ' In Pfaf_Over got INFO1 = ',info1,' from ZGETRF'
	stop
      end if
     
      !Inverse
      call ZGETRI(N,ZUbar_inv,N,IPIV,ZWORK,N,INFO2)
      if(info2.ne.0) then
	write(6,*) ' In Pfaf_Over got INFO2 = ',info2,' from ZGETRI'
	stop
      end if
      
 
      !Z
      call ZGEMM ('n','n',N,N,N,zone,ZVbar_c,N,ZUbar_inv,N,zzero,ZM1,N)
     
      !Big M
      do ii=1,N
       do jj=1,N
         ZPfaf(jj,ii)=ZM1(jj,ii)
         ZPfaf(jj+N,ii+N)=-ZM0(jj,ii)
         ZPfaf(jj,ii+N)=zzero
         ZPfaf(jj+N,ii)=zzero
         if (ii.eq.jj) then
         ZPfaf(ii,ii+N)=-zone
         ZPfaf(ii+N,ii)=zone
	 end if
       end do      
      end do
           
      call ZPfaffianF(ZPfaf,2*N,2*N,Ipiv2,ZdetPfaf)      
      
      end subroutine 












!+--------------------------------------------------------------------+
!|   C o m p l e x    V e r s i o n                                   |
!|                                                                    |
!+--------------------------------------------------------------------+
       Subroutine ZPfaffianF(SK,LDS,N,Ipiv,Pf)
       Implicit none
       Double Complex SK(LDS,N)
       Double Complex Pf,SS,SW,one,zero
       Double Precision epsln,phas,big
       Integer Ipiv(N,2),N,LDS,NB,IB,i,j,k,ip,NR,NC,I1,I2

       one = dcmplx(1.0d+00,0.0d+00)
       zero= dcmplx(0.0d+00,0.0d+00)
       
       epsln = 1.0d-13   ! smallest number such that 1+epsln=1
       
       if(mod(N,2).eq.1) then ! N odd
       
       Pf = zero
       
       else
          
       NB= N / 2   ! Number of 2x2 blocks
              
       if (NB.eq.1) then 

! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
       
       Pf = SK(1,2)
       
       else   ! NB.gt.1
       
       Pf = one
       
       do IB=1,NB-1
          
          NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
          NC = NR + 1

          big = 0.0d+00

          I1 = NR
          I2 = NC
                    
          do i=NR,N-1 ! all rows
             do j=i+1,N
          
                if(abs(SK( i , j )).gt. big) then
                   big = abs(SK( i , j ))
                   I1 = i
                   I2 = j 
                end if
                
             end do
          end do
          
          Ipiv( IB , 1) = I1 ! to initialize
          Ipiv( IB , 2) = I2 ! to initialize

          phas = 1.0d+00
          if(I1.eq.NR) phas = -phas
          if(I2.eq.NC) phas = -phas 
          
! Pivoting of element NR,NC with ip1,ip2

! 
!   Pivoting for a skew-symmetric matrix (Upper)
!
          if(I1.ne.NR) Call Zexch(SK,LDS,N,NR,I1)
          if(I2.ne.NC) Call Zexch(SK,LDS,N,NC,I2)
          
       ss = Sk( NR , NC )

       
!       write(6,*) ' IB I2, I2 ', IB,I1,I2,big,ss

       if(abs(ss).gt.epsln) then
!
!      Updating the Schur complement  matrix
!

       do i = 2*IB+1 , N-1
       
          do j = i+1 , N

             Sk(i,j) = Sk(i,j) +                               &
     &                (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
          end do
       end do
!
! Storing X and Y vectors in the lower part of the matrix
!       
       do i = 2*IB+1 , N
       
          Sk( i , NR ) =  -Sk( NC , i )/SS
          Sk( i , NC ) =   Sk( NR , i )/SS
          
       end do
        
       if (IB.ge.2) then  ! swap

!
!   Swapping
!       
       do j= 1 , 2*IB
       
          SW           = Sk( NR , j )
          Sk( NR , j ) = Sk( I1 , j ) 
          Sk( I1 , j ) = SW
       
          SW           = Sk( NC , j )
          Sk( NC , j ) = Sk( I2 , j ) 
          Sk( I2 , j ) = SW

       end do       
       
       end if
                      
       else
       
       Pf = zero
       
       return
       
       end if ! dabs(ss).gt.epsln
       
       Pf = Pf * SS * phas
       
       end do ! IB
       
       Pf = Pf * Sk(N-1,N)
       
       end if !    NB.eq.1
          
       end if !    mod(N,2).eq.1
       
       return
       end
!
!  Exchange of rows i1 i2 and columns i1 i2
!  SK is assumed upper triangular
!  i1 < i2
!
       Subroutine ZExch(SK,LDS,N,I1,I2)
       Implicit none
       Double Complex SK(LDS,N)
       Double Complex SS
       Integer LDS,N,I1,I2,I
              
       SK(i1,i2) = -SK(i1,i2)
                              
       if(I1.ne.1) then
          
          do i=1,I1-1
             SS           = SK( i , I1 )
             SK( i , I1 ) = SK( i , I2 )
             SK( i , I2 ) = SS
          end do
       end if

       if(I2.ne.N) then
          
          do i=I2+1,N
             SS           = SK( I1 , i )
             SK( I1 , i ) = SK( I2 , i )
             SK( I2 , i ) = SS
          end do
       end if

       if(I2.ge.I1+2) then
          
          do i=I1+1,I2-1
             SS           =  SK( I1 , i )
             SK( I1 , i ) = -SK( i , I2 )
             SK( i , I2 ) = -SS
          end do
       end if
       
       return
       end




!     ......................................................
        subroutine NZPROJ_PAV(U0,V0,akin,          &
     &  AZ_P,AN_P,                                 &
     &  EKin_P,Ekin_N,                             &
     &  EHF_PP,EHF_PN,EHF_NP,EHF_NN,               &
     &  EPa_PP,EPa_PN,EPa_NP,EPa_NN,               &
     &  EHFB_P,EHFB_N,NFOM,NLEV)
!     ......................................................

        USE VAPHFB_PAR
        implicit real*8 (a-h,o-y)
        implicit complex*16(z)

!        PARAMETER (NFOM=9)
!        PARAMETER (NOBS=15)

        DIMENSION U0(NLEV,NLEV)
        DIMENSION V0(NLEV,NLEV)
        DIMENSION ZU0(NLEV,NLEV)
        DIMENSION ZV0(NLEV,NLEV)
        DIMENSION ZU0bar(NLEV,NLEV)
        DIMENSION ZV0bar(NLEV,NLEV)
	DIMENSION ZU0tilde(NLEV,NLEV)
	DIMENSION ZV0tilde(NLEV,NLEV)
        DIMENSION ZZZc(NLEV,NLEV)
        DIMENSION ZZc(NLEV,NLEV)
        DIMENSION ZRO(NLEV,NLEV)
        DIMENSION Zkapa10(NLEV,NLEV)
        DIMENSION Zkapa01(NLEV,NLEV)
        DIMENSION zgamma(NLEV,NLEV)
        DIMENSION zdelta10(NLEV,NLEV)
	DIMENSION zdelta01(NLEV,NLEV)
        DIMENSION zham(NLEV,NLEV)
	


	DIMENSION ZAphi(NLEV,NLEV)
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DIMENSION akin(NLEV,NLEV)
        DIMENSION Zakin(NLEV,NLEV)
!        DIMENSION V_ABCD(NLEV,NLEV,NLEV,NLEV)
        
        Dimension Zsumphi(NOBS_PAV)
        
!        COMMON /particles/ Aprot,Aneut

        do ii=1,NOBS_PAV
         Zsumphi(ii)=zzero
        end do
        
        do ii=1,NLEV
         do jj=1,NLEV
          ZZc(jj,ii)   =zzero
          ZU0bar(jj,ii)=zzero
          ZV0bar(jj,ii)=zzero
          ZV0(jj,ii)=zone*V0(jj,ii)
          ZU0(jj,ii)=zone*U0(jj,ii)
          Zakin(jj,ii)=zone*akin(jj,ii)
         end do
        end do

!      commented by jmyao
      zDetvac=zone
      if(NFOM .gt. 1) then
        call zbareover(ZV0,ZU0,ZZc,Det,NLEV)
        zDetvac=dsqrt((1./Det))*zone
      endif 
	
	aaa=1.d0*NLEV*(NLEV+1)/2.d0
	snpfaf=(-1.d0)**(aaa)
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            FOMENKO STARTS                      !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO LLLP=0,NFOM-1   !Fomenko                
          
        phiLP=(2*pi/NFOM)*LLLP
        if(NFOM.eq.1) phiLP=0.d0*pi/180.d0

         zphas_p=cdexp(-Zimag*phiLP*Input%nprot)

        DO LLLN=0,NFOM-1   !Fomenko                
          
        phiLN=(2*pi/NFOM)*LLLN
        if(NFOM.eq.1) phiLN=0.d0*pi/180.d0
        
        zphas_n=cdexp(-Zimag*phiLN*Input%nneut)


	call UVphi(phiLP,phiLN,ZU0,ZV0,           &
     &	ZU0bar,ZV0bar,ZAphi,ZU0tilde,             &
     &	ZV0tilde,NLEV)                    

        
        !norm overlap

	call Pfaf_Over(ZZc,ZU0bar,ZV0bar,ZdetPfaf,NLEV)


       
!      ... jmyao
       znorm=snpfaf*zDetvac*ZdetPfaf*zphas_p*zphas_n

        if(NFOM.eq.1) znorm = zone
        zsumphi(1)=zsumphi(1)+znorm
       
!        write(*,*) 'PAV: znorm',znorm 

        !ro,kappa...
        call ROKAPPA_NZPROJ(ZU0,ZV0,ZU0tilde,ZV0tilde,      &
     &	ZRO,ZKAPA10,ZKAPA01,NLEV)	
        
        
        !fields...
        call HFB_FIELD_COMPLEX(zro,zkapa10,ZKAPA01,  &
     &  zakin,zgamma,zham,zdelta10,zdelta01,NLEV)  
    
       !energy
       call HFB_ENER_COMPLEX_fin(zro,zkapa01,zakin,   &
     &  zgamma,zdelta10,                              &
     &  ZEkin_P,ZEkin_N,                              &
     &  ZEHF_PP,ZEHF_PN,ZEHF_NP,ZEHF_NN,              &
     &  ZEPa_PP,ZEPa_PN,ZEPa_NP,ZEPa_NN,              &
     &  ZEHFB_P,ZEHFB_N,NLEV)


	zsumphi(2)=zsumphi(2)+ZEkin_P*znorm
	zsumphi(3)=zsumphi(3)+ZEkin_N*znorm
	zsumphi(4)=zsumphi(4)+ZEHF_PP*znorm
	zsumphi(5)=zsumphi(5)+ZEHF_PN*znorm
	zsumphi(6)=zsumphi(6)+ZEHF_NP*znorm
	zsumphi(7)=zsumphi(7)+ZEHF_NN*znorm
	zsumphi(8)=zsumphi(8)+ZEPa_PP*znorm
	zsumphi(9)=zsumphi(9)+ZEPa_PN*znorm
	zsumphi(10)=zsumphi(10)+ZEPa_NP*znorm
	zsumphi(11)=zsumphi(11)+ZEPa_NN*znorm
	zsumphi(12)=zsumphi(12)+ZEHFB_P*znorm
	zsumphi(13)=zsumphi(13)+ZEHFB_N*znorm
        

        !number of particles for testing purposes
	call ZN_0(ZRO,Zprot_0,Zneut_0,NLEV)
    	
        zsumphi(14)=zsumphi(14)+Zprot_0*znorm
        zsumphi(15)=zsumphi(15)+Zneut_0*znorm

      

        END DO !Fomenko
        END DO        



       EKin_P=dreal(zsumphi(2)/zsumphi(1)) 
       EKin_N=dreal(zsumphi(3)/zsumphi(1)) 
       EHF_PP=dreal(zsumphi(4)/zsumphi(1)) 
       EHF_PN=dreal(zsumphi(5)/zsumphi(1)) 
       EHF_NP=dreal(zsumphi(6)/zsumphi(1)) 
       EHF_NN=dreal(zsumphi(7)/zsumphi(1)) 
       EPa_PP=dreal(zsumphi(8)/zsumphi(1)) 
       EPa_PN=dreal(zsumphi(9)/zsumphi(1)) 
       EPa_NP=dreal(zsumphi(10)/zsumphi(1)) 
       EPa_NN=dreal(zsumphi(11)/zsumphi(1)) 
       EHFB_P=dreal(zsumphi(12)/zsumphi(1)) 
       EHFB_N=dreal(zsumphi(13)/zsumphi(1)) 


       AZ_P=dreal(zsumphi(14)/zsumphi(1))   
       AN_P=dreal(zsumphi(15)/zsumphi(1))   

!        print*,' ' 
!        write(6,'(A15,3F15.8)') 'norm',zsumphi(1)/(1.d0*NFOM**2)
!        write(6,'(A15,3F15.8)') 'Z',zsumphi(3)/zsumphi(1)
!        write(6,'(A15,3F15.8)') 'N',zsumphi(4)/zsumphi(1)
!        write(6,'(A15,3F15.8)') 'EPNP',zsumphi(2)/zsumphi(1)
        end subroutine









 
