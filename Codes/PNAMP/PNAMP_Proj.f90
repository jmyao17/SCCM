
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!! Norm ovelap with the bare vacuum !!!!!!!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

      Subroutine zbareover(ZV,ZU,ZZc,Det,N)
      Use VAPHFB_Par      
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
      
!      ulimit=1.e+15
      do ii=1,N
       do jj=1,N
        ZUINV(jj,ii)=ZU(jj,ii)
       end do
      end do


      !LU factorization: The factorization has the form A = P * L * U
      call ZGETRF(N,N,ZUINV,N,IPIV,INFO1)

      if(info1.ne.0) then
	write(6,*) ' In bareover got INFO1 = ',info1,' from ZGETRF'
	stop
      end if
    
!    This method inverts U and then computes inv(A) by solving the system
!    inv(A)*L = inv(U) for inv(A).
!    ............................. 
      ! U Inverse: U^-1
      call ZGETRI(N,ZUINV,N,IPIV,ZWORK,N,INFO2)
      if(info2.ne.0) then
	write(6,*) ' In bareover got INFO2 = ',info2,' from ZGETRI'
	stop
      end if
      

      ! Z^* = VU^-1
      call ZGEMM ('n','n',N,N,N,zone,ZV,N,ZUINV,N,zzero,ZZ,N)
     
      do ii=1,N      
      do jj=1,N
         ZZc(jj,ii)=dconjg(ZZ(jj,ii))
!         write(700,'(2i3,2f12.5)') ii,jj,ZZc(jj,ii)
      end do
      end do
      
      
      !Z^t*Z
      call ZGEMM ('t','n',N,N,N,zone,ZZc,N,ZZc,N,zzero,ZZTZ,N)
            
      !Eigensystem: ZGEEV computes for an N-by-N complex nonsymmetric matrix A, 
      !the eigenvalues and, optionally, the left and/or right eigenvectors.
      call ZGEEV('n','n',N,ZZTZ,N,ZWR,ZVL,N,ZVR,N,ZWORK1,2*N,WORK1,INFO3 )
          if(info3.ne.0) then
           write(6,*) ' In bareover got INFO3 = ',info3,' from ZGEEV'
             stop
      end if
      
      ztest = zzero
      zdet=zone
      do ii=1,N
       zdet=zdet*(zone+ZWR(ii))
       ztest = ztest + log(1.d0+ZWR(ii))
      end do
      Det  = dreal(EXP(-ztest/4.d0))
      return
      end  



      subroutine ROTATION_MAT(alp,bet,gam,ZROT,NLEV)
!    .................................................................
!     R(alp,bet,gam) = e^(-i alp J_z) * e^(-i bet J_y) * e^(-i gam J_z)
!    .................................................................
      USE VAPHFB_PAR
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      
      DIMENSION ZROT(NLEV,NLEV)
      DIMENSION ZAUX(NLEV,NLEV)
      do ii=1,NLEV
	do jj=1,NLEV
	 ZROT(jj,ii)=zzero
	end do
      end do
      
      !matrix elements
      do ia=1,NLEV

         ina=tnljm%n(ia) ! nlindex(ia)
         ila =tnljm%l(ia) ! lang(ia)
         ija =tnljm%twoj(ia) ! jang(ia)
         ima =tnljm%twom(ia) ! mjang(ia)
         ita =tnljm%t(ia)    ! mtisos(ia)

	rma=ima/2.d0

	zialpha_ma=cdexp(-Zimag*alp*rma) 

       do ib=1,NLEV

         inb =tnljm%n(ib) ! nlindex(ia)
         ilb =tnljm%l(ib) ! lang(ia)
         ijb =tnljm%twoj(ib) ! jang(ia)
         imb =tnljm%twom(ib) ! mjang(ia)
         itb =tnljm%t(ib)    ! mtisos(ia)

        if((inb.eq.ina).AND.(ilb.eq.ila).AND.  &
     &  (itb.eq.ita).AND.(ijb.eq.ija)) then

        rmb=imb/2.d0

        zigamma_mb=cdexp(-Zimag*gam*rmb)
        ZROT(ia,ib)=zialpha_ma*zigamma_mb*dwignerI_gen('h',ija,ima,imb,bet)
        end if 

       end do
      end do
     
!      call ZGEMM ('c','n',NLEV,NLEV,NLEV,zone, &     
!     & ZROT,NLEV,ZROT,NLEV,Zzero,ZAUX,NLEV)
!      call print_zmatrix(NLEV,ZROT)
      
      end subroutine





      SUBROUTINE Pfaf_Over(ZM0,ZUbar,ZVbar,ZdetPfaf,N)    
      USE VAPHFB_Par
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
!     .............. scaled by a factor 
      do ii=1,N*2
      do jj=1,N*2
           ZPfaf(ii,jj)=ZPfaf(ii,jj)*Input%Cscale
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

             Sk(i,j) = Sk(i,j) +   (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
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



      subroutine UVphiomega(philZ,philN,ZROT,ZU,ZV,ZUbar,ZVbar,ZAM_PHI_OMEGA,N)
!    ............................................................
!    calculate  U(Omega), V(Omega,phi_n,phi_p)
!    ............................................................

      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))
      
      Dimension ZV(N,N),ZU(N,N)
      Dimension ZUbar(N,N),ZVbar(N,N)
      Dimension ZROT(N,N)
      Dimension ZAM_PHI_OMEGA(N,N),ZAM_PHI_OMEGAc(N,N)

!     ............................
      zphasUZ=cdexp(zimag*philZ)       ! exp(i*phi_p)
      zphasUN=cdexp(zimag*philN)       ! exp(i*phi_n)
      zphasVZ=cdexp(-zimag*philZ)
      zphasVN=cdexp(-zimag*philN)
      
      do ii=1,N
       do jj=1,N
          ZAM_PHI_OMEGA(jj,ii)=zzero
          ZAM_PHI_OMEGAc(jj,ii)=zzero
        end do
      end do
      
!   .........................
!   U_p(Omega) = R(Omega)   U_p: 1,.., N/2
!   U_n(Omega) = R(Omega)   U_n: N/2+1, ..., N
!   V_p(Omega) = R^*(Omega) V_p: 1,..., N/2
!   V_n(Omega) = R^*(Omega) V_n: N/2+1, ..., N
!   ........................

      do ii=1,N/2
       do jj=1,N/2
        ZAM_PHI_OMEGA(ii,jj)         =zphasUZ*ZROT(ii,jj)
        ZAM_PHI_OMEGA(ii+N/2,jj+N/2) =zphasUN*ZROT(ii+N/2,jj+N/2)
        ZAM_PHI_OMEGAc(ii,jj)        =dconjg(ZAM_PHI_OMEGA(ii,jj))
        ZAM_PHI_OMEGAc(ii+N/2,jj+N/2)=dconjg(ZAM_PHI_OMEGA(ii+N/2,jj+N/2))
       end do
      end do

      !Ubar
      call ZGEMM ('n','n',N,N,N,zone,ZAM_PHI_OMEGA,N,ZU,N,zzero,ZUbar,N)
!    ............................................
!     If you find ZUbar  different, which is probabbly due to different ZROT
!     from the different definition of the ordering for the s.p. energy levels
!    ............................................
!      write(*,*) N
!      call ZGEMM ('n','n',N,N,N,zone,ZROT,N,ZU,N,zzero,ZUbar,N)
!      call print_zmatrix3(N,ZROT,ZU,ZUbar)  ! ZROT, ZU are the same

      !Vbar
      call ZGEMM ('n','n',N,N,N,zone,ZAM_PHI_OMEGAc,N,ZV,N,zzero,ZVbar,N)

           
      end subroutine



       SUBROUTINE ROKAPPAphiomega(ZAM_phi_omega,ZU0,ZV0,ZU1,ZV1,ZRO,ZKAPA10,ZKAPA01,N)
     
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zonem=(-1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))

      DIMENSION ZAM_phi_omega(N,N)
      DIMENSION ZAM_phi_omegac(N,N)
      DIMENSION ZU0(N,N),ZV0(N,N)
      DIMENSION ZU1(N,N),ZV1(N,N)
      DIMENSION ZRO(N,N)
      DIMENSION ZKAPA10(N,N),ZKAPA01(N,N)

      DIMENSION ZT22INV(N,N)
      Dimension IPIV(N)
      Dimension ZWORK(N)

      DIMENSION ZAUX1(N,N)
      DIMENSION ZAUX2(N,N)
      DIMENSION ZAUX3(N,N)
      DIMENSION ZAUX4(N,N)


     
      
      do ii=1,N
       do jj=1,N
	ZAM_phi_omegac(jj,ii)=dconjg(ZAM_phi_omega(jj,ii))
        end do
      end do
     
!!!!!!!!!!!!! T22/T22inv !!!!!!!!!!!!!!!!!     
     
      
      call ZGEMM ('t','n',N,N,N,zone,        &
     &	ZV0,N,ZAM_phi_omega,N,zzero,ZAUX1,N)

      call ZGEMM ('n','n',N,N,N,zone,        &
     &	ZAUX1,N,ZV1,N,zzero,ZAUX2,N)
     
       call ZGEMM ('t','n',N,N,N,zone,       &
     &	ZU0,N,ZAM_phi_omegac,N,zzero,ZAUX3,N)

      call ZGEMM ('n','n',N,N,N,zone,        &
     &	ZAUX3,N,ZU1,N,zzero,ZAUX4,N)
    
     
      do ii=1,N
       do jj=1,N
          ZT22INV(jj,ii)=ZAUX2(jj,ii)+ZAUX4(jj,ii)
        end do
      end do
     
      !LU factorization
      call ZGETRF(N,N,ZT22INV,N,IPIV,INFO1)
      !Inverse
      call ZGETRI(N,ZT22INV,N,IPIV,ZWORK,N,INFO2)
     

!!!!!!!!!!! RHO, KAPPA !!!!!!!!!!!
      call ZGEMM ('n','n',N,N,N,zone,  &
     &	ZAM_phi_omega,N,ZV1,N,zzero,ZAUX1,N)
     

      call ZGEMM ('n','n',N,N,N,zonem,  &
     &	ZAM_phi_omegac,N,ZU1,N,zzero,ZAUX2,N)
     

      call ZGEMM ('n','t',N,N,N,zone,   &
     &	ZT22INV,N,ZV0,N,zzero,ZAUX3,N)

      call ZGEMM ('n','t',N,N,N,zone,   &
     &	ZT22INV,N,ZU0,N,zzero,ZAUX4,N)

      !rho
      call ZGEMM ('n','n',N,N,N,zone,  &
     &	ZAUX1,N,ZAUX3,N,zzero,ZRO,N)

      !kappa10
      call ZGEMM ('n','n',N,N,N,zone,   &
     &	ZAUX1,N,ZAUX4,N,zzero,ZKAPA10,N)

      !kappa10
      call ZGEMM ('n','n',N,N,N,zone,  &
     &	ZAUX2,N,ZAUX3,N,zzero,ZKAPA01,N)


     
       return
       END SUBROUTINE 


      subroutine UVbar(philZ,philN,ZU,ZV,ZUbar,ZVbar,N)
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))
      
      Dimension ZV(N,N),ZU(N,N)
      Dimension ZUbar(N,N),ZVbar(N,N)
      
      zphasUZ=cdexp(zimag*philZ)
      zphasUN=cdexp(zimag*philN)
      zphasVZ=cdexp(-zimag*philZ)
      zphasVN=cdexp(-zimag*philN)
      
      do ii=1,N/2
       do jj=1,N
        ZUbar(ii    ,jj)=zphasUZ*ZU(ii    ,jj)
        ZUbar(ii+N/2,jj)=zphasUN*ZU(ii+N/2,jj)
        ZVbar(ii    ,jj)=zphasVZ*ZV(ii    ,jj)
        ZVbar(ii+N/2,jj)=zphasVN*ZV(ii+N/2,jj)
       end do
      end do
     
      
      
      end subroutine


      subroutine UVphi(ZU,ZV,ZUbar,ZVbar,ZAphi,ZRO,ZKAPA10,ZKAPA01,ZUphic,ZVphic,N)

      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      Parameter (zone=(1.d0,0.d0))
      Parameter (zonem=(-1.d0,0.d0))
      Parameter (zimag=(0.d0,1.d0))
      Parameter (zzero=(0.d0,0.d0))
      
      Dimension ZV(N,N),ZU(N,N)
      Dimension ZVc(N,N),ZUc(N,N)
      Dimension ZUbar(N,N),ZVbar(N,N)
      Dimension ZUphi(N,N),ZVphi(N,N)
      Dimension ZUphic(N,N),ZVphic(N,N)
      Dimension ZUtilde(N,N),ZVtilde(N,N)
      Dimension ZUtilde_inv(N,N),ZVtilde_c(N,N)
      Dimension ZAphi(N,N)
      Dimension ZRO(N,N),ZKAPA10(N,N),ZKAPA01(N,N)
      Dimension IPIV(N)
      Dimension ZWORK(N)
      Dimension ZAUX1(N,N)
      Dimension ZAUX2(N,N)

      call ZGEMM ('c','n',N,N,N,zone,ZU,N,ZUbar,N,Zzero,ZUtilde,N)	 
      call ZGEMM ('c','n',N,N,N,zone,ZV,N,ZVbar,N,zone,ZUtilde,N)	 

      call ZGEMM ('t','n',N,N,N,zone,ZV,N,ZUbar,N,Zzero,ZVtilde,N)	 
      call ZGEMM ('t','n',N,N,N,zone,ZU,N,ZVbar,N,zone,ZVtilde,N)
     
      do ii=1,N
       do jj=1,N
        ZUtilde_inv(jj,ii)=dconjg(ZUtilde(jj,ii))
        ZVtilde_c(jj,ii)=dconjg(ZVtilde(jj,ii))
        ZVc(jj,ii)=dconjg(ZV(jj,ii))
        ZUc(jj,ii)=dconjg(ZU(jj,ii))
       end do
      end do	 

      Call ZGETRF(N,N,ZUtilde_inv,N,IPIV,INFO)

      if(info.ne.0) then
        write(6,*) ' In AINV got INFO = ',info,' from ZGETRF 3'
        stop
      end if
      
      Call Zgetri (N,ZUtilde_inv,N,IPIV,ZWORK,N,INFO) 
      if(info.ne.0) then
        write(6,*) ' In AINV got INFO = ',info,' from DGETRI 4'
        stop
      end if

      call ZGEMM ('n','n',N,N,N,Zone,          &
     &	ZVtilde_c,N,ZUtilde_inv,N,zzero,ZAphi,N)
     
      call ZGEMM ('n','n',N,N,N,zone,	      &
     & ZU,N,ZAphi,N,Zzero,ZAUX1,N)	 
      call ZGEMM ('n','n',N,N,N,zone,	      &
     & ZV,N,ZAphi,N,Zzero,ZAUX2,N)
     
     
      do ii=1,N
       do jj=1,N
        ZVphi(jj,ii)=ZVc(jj,ii)+ZAUX1(jj,ii)
        ZUphi(jj,ii)=ZUc(jj,ii)+ZAUX2(jj,ii)
        ZVphic(jj,ii)=dconjg(ZVphi(jj,ii))
        ZUphic(jj,ii)=dconjg(ZUphi(jj,ii))
       end do
      end do	 
      
      
      call ZGEMM ('n','t',N,N,N,zone, ZVphi,N,ZV,N,Zzero,ZRO,N)
      
      call ZGEMM ('n','t',N,N,N,zone, ZVphi,N,ZU,N,Zzero,ZKAPA10,N)

      call ZGEMM ('n','t',N,N,N,zonem, ZUphi,N,ZV,N,Zzero,ZKAPA01,N)
      
      end subroutine


      
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Computes the norm overlap following Neergard's prescription           c
!                                                                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Neer_Over(Zneer,ZUbar,ZVbar,ZdetNe,N)    

      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
      Parameter (Zone=(1.d0,0.d0))
      Parameter (Zzero=(0.d0,0.d0))
      Parameter (Zimag=(0.d0,1.d0))

      Dimension Zneer(N,N) !rotation matrix 
      Dimension ZUbar(N,N),ZVbar(N,N) 
      Dimension ZUbar_inv(N,N),ZVbar_c(N,N) 
      Dimension ZZ(N,N)
      
      Dimension IPIV(N)
      Dimension ZWORK(N)

      Dimension ZZZS(N,N) 
     
      Dimension ZAUX1(N,N) 
      Dimension ZAUX2(N,N) 
      
      do ii=1,N
       do jj=1,N
        ZUbar_inv(jj,ii)=dconjg(ZUbar(jj,ii))
	ZVbar_c(jj,ii)=dconjg(ZVbar(jj,ii))
       end do
      end do
      !LU factorization
      call ZGETRF(N,N,ZUbar_inv,N,IPIV,INFO1)


      if(info1.ne.0) then
	write(6,*) ' In Neerover got INFO1 = ',info1,' from ZGETRF'
	stop
      end if
     
      !Inverse
      call ZGETRI(N,ZUbar_inv,N,IPIV,ZWORK,N,INFO2)
      if(info2.ne.0) then
	write(6,*) ' In Neerover got INFO2 = ',info2,' from ZGETRI'
	stop
      end if
      

      !Z
      call ZGEMM ('n','n',N,N,N,zone,ZVbar_c,N,ZUbar_inv,N,zzero,ZZ,N)
           
      
      !ZtZ
      call ZGEMM ('t','n',N,N,N,zone,Zneer,N,ZZ,N,zzero,ZZZS,N)
     
            

      
      call ccstot(ZZZS,ZdetNe,N) 
       
      return
      
      END SUBROUTINE     
 
 
      SUBROUTINE ccstot(ZCCS,ZdetNe,N)
      
      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
      Parameter (Zone=(1.d0,0.d0))
      Parameter (Zzero=(0.d0,0.d0))
      Parameter (Zimag=(0.d0,1.d0))
      
      Dimension ZCCS(N,N) 
      Dimension Zeigen(N),Zeigenaux1(N)
      Dimension Zeigenaux(N)
      Dimension ZWORK(2*N)
      Dimension AModeigen(N)
      Dimension Ieigen(N)
      
      Dimension RWORK(2*N)
      
      Dimension ZVL(N,N),ZVR(N,N)

!      !Diagonalization
      Call ZGEEV('n','n',N,ZCCS,N,Zeigen,ZVL,N,ZVR,N,   &
     & ZWORK,2*N,RWORK,INFOCC)
      if(infoCC.ne.0) then
       write(6,*) ' In ccstot got INFOCC = ',infoCC,' from ZGETRF CC'
       stop
      end if
 
      do ii=1,N
       Amodeigen(ii)=dsqrt(dreal(Zeigen(ii))**2+dimag(Zeigen(ii))**2)
       Ieigen(ii)=ii
      end do
      
      call SORTVC (N,Ieigen,Amodeigen,1)
      
      
      do ii=1,N
       jj=Ieigen(ii)
       Zeigenaux1(ii)=Zeigen(jj)
      end do
      
      do ii=1,N,2

       Zeigenaux(ii)=Zeigenaux1(ii)
       Zeigenaux(ii+1)=Zeigenaux1(ii+1)
       
       aa=Dimag(Zeigenaux1(ii));bb=Dimag(Zeigenaux1(ii+1))
       if(aa.lt.bb) then 
        Zeigenaux(ii+1)=Zeigenaux1(ii)
        Zeigenaux(ii)=Zeigenaux1(ii+1)
       end if
      end do
      
      ZdetNe=Zone
      ij=0
      do ii=1,N,2
       jj=ii+ij
       ZdetNe=ZdetNe*(Zone+Zeigenaux(jj)) 
       ij=1-ij  
      end do
      


1000  format(I4,2E20.10)
      return
      
      END SUBROUTINE      
      
! ******************************************************************
!
!  Sortvc sorts D and (if required) E and the columns of Q.
!
!     Prameters:
!
!       ND      (I) :  Number of elements in D (or columns in Q)
!       ID      (IO):  original index corresponding to the elements sorted 
!       D       (I) :  Vector to sort
!               (O) :  Sorted vector        
!       M       (I) :  1: Descending order in D
!                     -1: Ascending order in D
!                      otherwise: no sorting
! **********************************************************************
!
      SUBROUTINE SORTVC (ND,ID,D,M)
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL    LMIN,LMAX
        DIMENSION  D(ND)
	DIMENSION  ID(ND)
!
        IF (ND .LT. 2) RETURN
        LMAX = (M .EQ.  1)
        LMIN = (M .EQ. -1)
        IF (.NOT. (LMAX .OR. LMIN)) RETURN
        DO 40 KK = 2,ND
          K = KK - 1
          J = K
          H = D(K)
	  L = ID(K)
!
!  FIND EXTREMUM
!
          DO 10 I = KK,ND
            S = D(I)
	    MM = ID(I)
            IF (LMIN .AND. (S .GE. H)) GOTO 10
            IF (LMAX .AND. (S .LE. H)) GOTO 10
            J = I
            H = S
	    L = MM
   10     CONTINUE
          IF (J .EQ. K) GOTO 40
!
!  SORT D
!
          D(J) = D(K)
	  ID(J) = ID(K)
          D(K) = H
	  ID(K) = L

   40   CONTINUE
        RETURN
      END
           
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!! ISOSPIN ROTATION MATRIX  !!!!!!!!!!!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine isospin_rot(alp,bet,gam,ZROT_ISO,NDIM)

	implicit real*8 (a-h,o-y)
	implicit complex*16 (z)
	PARAMETER (zone=(1.d0,0.d0))
	PARAMETER (zzero=(0.d0,0.d0))
	PARAMETER (zimag=(0.d0,1.d0))
	
	DIMENSION ZROT_ISO(NDIM,NDIM)
	
	
	do ii=1,NDIM
	 do jj=1,NDIM
	  ZROT_ISO(jj,ii)=zero
	 end do
	end do
	
	zprot_a=cdexp(zimag*alp/2.d0)
	zprot_c=cdexp(zimag*gam/2.d0)

	zneut_a=cdexp(-zimag*alp/2.d0)
	zneut_c=cdexp(-zimag*gam/2.d0)
	
	sb=dsin(bet/2.)
	cb=dcos(bet/2.)
	
!POINTERS
        ip=0
	ipb=ip+NDIM/4
	in=ipb+NDIM/4
	inb=in+NDIM/4   
	
	do ii=1,NDIM/4
	  !Tz
	  ZROT_ISO(ii+ip,ii+ip )  = zprot_a*zprot_c*cb
	  ZROT_ISO(ii+ipb,ii+ipb )= zprot_a*zprot_c*cb
	  ZROT_ISO(ii+in,ii+in )  = zneut_a*zneut_c*cb
	  ZROT_ISO(ii+inb,ii+inb )= zneut_a*zneut_c*cb
     
	  !Tx
	  ZROT_ISO(ii+ip,ii+in)  = zprot_a*zneut_c*sb
	  ZROT_ISO(ii+ipb,ii+inb)= zprot_a*zneut_c*sb
	  ZROT_ISO(ii+in,ii+ip)  =-zneut_a*zprot_c*sb
	  ZROT_ISO(ii+inb,ii+ipb)=-zneut_a*zprot_c*sb	 
	end do     
        
     
        end subroutine
	
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!! SPIN ROTATION MATRIX     !!!!!!!!!!!		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine spin_rot(alp,bet,gam,ZROT_SPIN,NDIM)

	implicit real*8 (a-h,o-y)
	implicit complex*16 (z)
	PARAMETER (zone=(1.d0,0.d0))
	PARAMETER (zzero=(0.d0,0.d0))
	PARAMETER (zimag=(0.d0,1.d0))
	
	DIMENSION ZROT_SPIN(NDIM,NDIM)
	
	
	do ii=1,NDIM
	 do jj=1,NDIM
	  ZROT_SPIN(jj,ii)=zero
	 end do
	end do
	
	zalpp=cdexp(zimag*alp/2.d0)
	zgamp=cdexp(zimag*gam/2.d0)
	zalpm=cdexp(-zimag*alp/2.d0)
	zgamm=cdexp(-zimag*gam/2.d0)
	
	sb=dsin(bet/2.)
	cb=dcos(bet/2.)
	
!POINTERS
        ip=0
	ipb=ip+NDIM/4
	in=ipb+NDIM/4
	inb=in+NDIM/4   
	
	do ii=1,NDIM/4
	  !protons
	  ZROT_SPIN(ii+ip,ii+ip )  = zalpm*cb*zgamm
	  ZROT_SPIN(ii+ipb,ii+ipb )= zalpp*cb*zgamp
	  ZROT_SPIN(ii+in,ii+in )  = zalpm*cb*zgamm
	  ZROT_SPIN(ii+inb,ii+inb )= zalpp*cb*zgamp
     
	  !neutrons
	  ZROT_SPIN(ii+ip,ii+ipb)=-zalpm*sb*zgamp
	  ZROT_SPIN(ii+ipb,ii+ip)=zalpp*sb*zgamm
	  ZROT_SPIN(ii+in,ii+inb)=-zalpm*sb*zgamp
	  ZROT_SPIN(ii+inb,ii+in)=zalpp*sb*zgamm    
	end do     
        
     
        end subroutine      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Evaluates:                                             c
!                                                            c
!     a) T22                                                 c
!                                                            c
!     For triaxial PNAMP                                     c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      
       SUBROUTINE T22_Over(phil,ZU1,ZV1,ZROT,ZT22_INV,zdeterm,N)


      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
      Parameter (Zone=(1.d0,0.d0))
      Parameter (Zonem=(-1.d0,0.d0))
      Parameter (Zzero=(0.d0,0.d0))
      Parameter (Zimag=(0.d0,1.d0))

      Dimension ZU1(N,N)
      Dimension ZV1(N,N)    
      Dimension ZU1c(N,N)
      Dimension ZV1c(N,N)    
      Dimension ZROT(N,N)
      Dimension ZROTconj(N,N)
      Dimension ZT22_INV(N,N)

      Dimension IPIV(N)
      Dimension ZWORK(N)

      
      Dimension ZAUX1(N,N)
      Dimension ZAUX2(N,N)
      
      Zfas1=cdexp(Zimag*phil)
      Zfas2=cdexp(-Zimag*phil)
      pi=4.d0*datan(1.d0)
      
      do ii=1,N
       do jj=1,N
        ZROTconj(jj,ii)=dconjg(ZROT(jj,ii))
        ZU1c(jj,ii)=dconjg(ZU1(jj,ii))
        ZV1c(jj,ii)=dconjg(ZV1(jj,ii))
       end do
      end do
       
      call ZGEMM ('t','n',N,N,N,Zfas1,ZV1,N,ZROT,N,Zzero,ZAUX1,N)
      call ZGEMM ('n','n',N,N,N,Zone,ZAUX1,N,ZV1c,N,Zzero,ZT22_INV,N)
   
      call ZGEMM ('t','n',N,N,N,Zfas2,ZU1,N,ZROTconj,N,Zzero,ZAUX2,N)
      call ZGEMM ('n','n',N,N,N,Zone,ZAUX2,N,ZU1c,N,Zone,ZT22_INV,N)
     
     
           
      Call ZGETRF(N,N,ZT22_INV,N,IPIV,INFO)

      if(info.ne.0) then
        write(6,*) ' In AINV got INFO = ',info,' from ZGETRF 3'
        stop
      end if
      
      zdeterm=zone
      inter=0
      do ii=1,N
       zdeterm=zdeterm*ZT22_INV(ii,ii)
         if(ipiv(ii).ne.ii) inter = inter + 1
      end do
      
      isign=1-2*Mod(inter,2)
      If(isign.eq.1) then
       zdeterm = zdeterm
      else
        zdeterm = -zdeterm
      end if
      
      zdeterm=zdeterm**(.5d0)
      


      Call Zgetri (N,ZT22_INV,N,IPIV,ZWORK,N,INFO) 
      if(info.ne.0) then
        write(6,*) ' In AINV got INFO = ',info,' from DGETRI 4'
        stop
      end if

      return
      END SUBROUTINE
      

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Computes the norm overlap following Neergard's prescription           c
!                                                                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Neer_Over_ISOS(phil,ZTHOUL,ZROT,ZdetNe,N)    

      implicit real*8(a-h,o-y)
      implicit complex*16 (z)
      
      Parameter (Zone=(1.d0,0.d0))
      Parameter (Zzero=(0.d0,0.d0))
      Parameter (Zimag=(0.d0,1.d0))

      Dimension ZROT(N,N) !rotation matrix 
      Dimension ZTHOUL(N,N) 
      
      Dimension ZZZS(N,N) 
     
      Dimension ZAUX1(N,N) 
      Dimension ZAUX2(N,N) 
      
      Zfas=cdexp(2.d0*Zimag*phiL)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!                  Calculation of auxiliary matrices                                        c
!                     									    c   
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call ZGEMM ('n','t',N,N,N,Zfas, ZTHOUL,N,ZROT,N,Zzero,ZAUX1,N)	 

      call ZGEMM ('t','n',N,N,N,Zone,ZTHOUL,N,ZROT,N,Zzero,ZAUX2,N)	 
      
      call ZGEMM ('n','n',N,N,N,Zone, ZAUX1,N,ZAUX2,N,Zzero,ZZZS,N)	



      
      call ccstot(ZZZS,ZdetNe,N) 
       
      return
      
      END SUBROUTINE    
