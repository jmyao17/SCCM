
        subroutine zMatrix_Rho3B(k1,k2,k3,k4,k5,k6,ZRO,ZKAPA10,ZKAPA01,ZRho3B,NLEV)
        USE VAPHFB_PAR
        implicit none
        integer NLEV
        integer k1,k2,k3,k4,k5,k6
        complex*16 ZRO,ZKAPA10,ZKAPA01,zMatrix,ZRho3B
        DIMENSION ZRO(NLEV,NLEV)    ! density matrix
        DIMENSION ZKAPA10(NLEV,NLEV)    ! density matrix
        DIMENSION ZKAPA01(NLEV,NLEV)    ! density matrix
        dimension zMatrix(6,6)
        Integer Ipiv2(6,2)
!       ................. initialization
        zMatrix =zzero 
!       ................. check
        if(max(k1,k2,k3,k4,k5,k6).gt.NLEV) &
     &  stop 'k values are not correct'

        zMatrix(1,2) = zkapa01(k1,k2)
        zMatrix(1,3) = zkapa01(k1,k3)
        zMatrix(2,3) = zkapa01(k2,k3)

        zMatrix(1,4) = zro(k6,k1)
        zMatrix(1,5) = zro(k5,k1)
        zMatrix(1,6) = zro(k4,k1)

        zMatrix(2,4) = zro(k6,k2)
        zMatrix(2,5) = zro(k5,k2)
        zMatrix(2,6) = zro(k4,k2)

        zMatrix(3,4) = zro(k6,k3)
        zMatrix(3,5) = zro(k5,k3)
        zMatrix(3,6) = zro(k4,k3)

        zMatrix(4,5) = zkapa10(k5,k6)
        zMatrix(4,6) = zkapa10(k4,k6)
        zMatrix(5,6) = zkapa10(k4,k5)

        call ZPfaffianF(zMatrix,6,6,Ipiv2,ZRho3B)      

        end  

!     .....................................................................      
!     Computing the three-body density in J-scheme
!     .....................................................................      
        subroutine ZRho3B_Angles(iangle,ZRO,ZKAPA10,ZKAPA01,znorm,ZRO3BJJ,NLEV)
        USE VAPHFB_PAR
        implicit none
        integer   NLEV,iangle 
        complex*16 ZRO,ZKAPA10,ZKAPA01,Znorm,Znorm2,ZRho3B
        DIMENSION ZRO(NLEV,NLEV)        ! density matrix           <c+c>
        DIMENSION ZKAPA10(NLEV,NLEV)    ! ab-normal density matrix <cc>
        DIMENSION ZKAPA01(NLEV,NLEV)    ! ab-normal density matrix <c+c+>
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRO3BJJ  

        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER k1,t1,n1,lj1,l1,j1,twom1 
        INTEGER k2,t2,n2,lj2,l2,j2,twom2 
        INTEGER k3,t3,n3,lj3,l3,j3,twom3 
        INTEGER k4,t4,n4,lj4,l4,j4,twom4 
        INTEGER k5,t5,n5,lj5,l5,j5,twom5 
        INTEGER k6,t6,n6,lj6,l6,j6,twom6
        INTEGER J123,P123,P456,J12,J45,P12,P45,t12,t45
        INTEGER id123,id456,b12,a12,b45,a45,M12,M45,M123 
        INTEGER tnlj1,tnlj2,tnlj4,tnlj5
        INTEGER SPB_l_lj,SPB_twoj_lj
        INTEGER icheck,id456_m,iUse_Sym

        REAL(DP) CG
        real*8 cg1,cg2,cg3,cg4,cg21,cg41
        real*8 omp_get_wtime,time0,time1
        complex*16 Trace
        real*8  Metric

        iUse_Sym=0
!     .............        
!        if(.NOT. ALLOCATED(ZRO3BJJ))  &
!     &  ALLOCATE(ZRO3BJJ(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax))

        Trace =zzero
        Znorm2 = Znorm/(PNP%NFOM**2)
        time0 = omp_get_wtime()
        write(*,'(a30)',advance='no') '  Computing the Rho3B ....'
!       write(*,*) Rho3B%idx123_vmax

        icheck = Rho3B%idx123_vmax/10


!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Trace,icheck,zkapa01,zro,zkapa10,&
!$OMP                                  & tnljm,CG_Save,Input,iv,Znorm2,Rho3B,ZRO3BJJ) 
!$OMP DO schedule(dynamic)  REDUCTION(+:Trace) 
        do id123=1,Rho3B%idx123_vmax

!          if(mod(id123,icheck).eq.0) call progress(id123/icheck)
           J123 = Rho3B%J123(id123)            ! doubled

           J12 = Rho3B%J12(id123)
           P12 = Rho3B%P12(id123)
           t12 = Rho3B%t12(id123)
           t3  = Rho3B%t3(id123)  
           lj3 = Rho3B%lj3(id123)  
           n3  = Rho3B%n3(id123)  
           l3  = SPB_l_lj(lj3)*2     ! doubled l3
           j3  = SPB_twoj_lj(lj3)    ! doubled
           P123 = mod(l3/2+P12,2)

           b12   = Rho3B%VTPB%block(J12,P12,t12)
           a12   = Rho3B%a12(id123) 
!      ............. quantum numbers for 1 and 2
           n1  = Rho3B%VTPB%n1(b12,a12) 
           t1  = Rho3B%VTPB%t1(b12,a12) 
           L1  = Rho3B%VTPB%twol1(b12,a12)  ! doubled 
           j1  = Rho3B%VTPB%j1(b12,a12)  ! doubled
           lj1 = (L1+j1-1)/2
           n2  = Rho3B%VTPB%n2(b12,a12)
           t2  = Rho3B%VTPB%t2(b12,a12)
           L2  = Rho3B%VTPB%twol2(b12,a12)  ! doubled 
           j2  = Rho3B%VTPB%j2(b12,a12)  ! doubled
           lj2 = (L2+j2-1)/2
!      ............. check
           tnlj1 = Rho3B%VSPB%tnlj(lj1,n1,t1)
           tnlj2 = Rho3B%VSPB%tnlj(lj2,n2,t2)
           if(a12 .ne. Rho3B%VTPB%a(b12,tnlj1,tnlj2)) &
     &     stop 'Error in QNs for a(1,2)'


!      ............... symmetry
        do id456=1,Rho3B%idx123_vmax

!        if(iUse_Sym.eq.1) id456_m = id123 
!        do id456=id456_m, Rho3B%idx123_vmax

           J45 = Rho3B%J12(id456)
           P45 = Rho3B%P12(id456)
           t45 = Rho3B%t12(id456)
           t6  = Rho3B%t3(id456)  
           lj6 = Rho3B%lj3(id456)
           n6  = Rho3B%n3(id456) 
           l6  = SPB_l_lj(lj6)*2       ! doubled l6
           j6  = SPB_twoj_lj(lj6)
           P456 = mod(l6/2+P45,2)

!      ........................................ check this symmetry
           if(t12+t3.ne.t45+t6) cycle 
           if(P123.ne.P456)     cycle 
           if(J123.ne.Rho3B%J123(id456)) cycle 

           b45   = Rho3B%VTPB%block(J45,P45,t45)
           a45   = Rho3B%a12(id456) 
!      ............. quantum numbers for 1 and 2
           n4  = Rho3B%VTPB%n1(b45,a45)
           t4  = Rho3B%VTPB%t1(b45,a45)
           L4  = Rho3B%VTPB%twol1(b45,a45)  ! doubled 
           j4  = Rho3B%VTPB%j1(b45,a45)     ! doubled
           lj4 = (L4+j4-1)/2
           n5  = Rho3B%VTPB%n2(b45,a45)
           t5  = Rho3B%VTPB%t2(b45,a45)
           L5  = Rho3B%VTPB%twol2(b45,a45)  ! doubled 
           j5  = Rho3B%VTPB%j2(b45,a45)  ! doubled
           lj5 = (L5+j5-1)/2
!      ............. check
!           tnlj4 = Rho3B%VSPB%tnlj(lj4,n4,t4)
!           tnlj5 = Rho3B%VSPB%tnlj(lj5,n5,t5)
!           if(a45 .ne. Rho3B%VTPB%a(b45,tnlj4,tnlj5)) &
!     &     stop 'Error in QNs for aa(4,5)'

          do twom1=-j1,j1,2
          do twom2=-j2,j2,2
             M12  = (twom1+twom2)/2
             cg1  = CG_Save((j1+1)/2,(twom1+1)/2,(j2+1)/2,(twom2+1)/2,J12,M12) ! J12,M12 are NOT doubled 

             if(abs(cg1).lt.CHOP) cycle

          do twom3=-j3,j3,2
             M123 = twom1+twom2+twom3 

!       ................................................
!       caution: 
!       One cannot use the following symmetry because the projection
!       operator is only commutating with the whole density
!       operator, not J-component.
!       ................................................
!             if(M123.ne.J123) cycle     ! only pick up the component with M=J
!       .............................................
             cg2 =  iv((j3+twom3)/2)*dsqrt((J123+1.d0)/(2*J12+1.d0)) &
     &       *CG_Save((J123+1)/2,(-M123+1)/2,(j3+1)/2,(twom3+1)/2,J12,-M12)  

!             cg21  = CG(J12*2,M12*2,j3,twom3,J123,M123)
!            if(abs(cg2-cg21).gt.CHOP) stop 'Error: dim. in CG_Save is &
!                                         &  NOT sufficient ! '

             if(Input%IntType.eq.0 .and. (J123+1)/2 .gt. 14) then
               write(*,*) 'PNAMP_Rho3B: 14 is smaller than (J123+1)/2',J123
               stop 
             endif

             if(abs(cg2).lt.CHOP) cycle


          do twom4=-j4,j4,2
          do twom5=-j5,j5,2
             M45  = (twom4+twom5)/2
             

             cg3 = CG_Save((j4+1)/2,(twom4+1)/2,(j5+1)/2,(twom5+1)/2,J45,M45) ! J12,M12 are NOT doubled 

             if(abs(cg3).lt.CHOP) cycle

             twom6 = M123-2*M45

             cg4 =  iv((j6+twom6)/2)*dsqrt((J123+1.d0)/(2*J45+1.d0)) &
     &              *CG_Save((J123+1)/2,(-M123+1)/2,(j6+1)/2,(twom6+1)/2,J45,-M45)  

             if(abs(cg4).lt.CHOP) cycle

!     ................
           k1 = tnljm%level(twom1,lj1,n1+1,t1)  ! in tnljm%level, n starts from 1
           k2 = tnljm%level(twom2,lj2,n2+1,t2)  ! in tnljm%level, n starts from 1
           k3 = tnljm%level(twom3,lj3,n3+1,t3)  ! in tnljm%level, n starts from 1
           k4 = tnljm%level(twom4,lj4,n4+1,t4)  ! in tnljm%level, n starts from 1
           k5 = tnljm%level(twom5,lj5,n5+1,t5)  ! in tnljm%level, n starts from 1
           k6 = tnljm%level(twom6,lj6,n6+1,t6)  ! in tnljm%level, n starts from 1


!          call zMatrix_Rho3B(k1,k2,k3,k4,k5,k6,ZRO,ZKAPA10,ZKAPA01,ZRho3B)

           ZRho3B = zkapa01(k1,k2)*(zro(k4,k3)*zkapa10(k6,k5) &
     &                             -zro(k5,k3)*zkapa10(k6,k4) &
     &                             +zro(k6,k3)*zkapa10(k5,k4))&
     &             -zkapa01(k1,k3)*(zro(k4,k2)*zkapa10(k6,k5) &
     &                             -zro(k5,k2)*zkapa10(k6,k4) &
     &                             +zro(k6,k2)*zkapa10(k5,k4))&
     &             +zro(k4,k1)*(zkapa01(k2,k3)*zkapa10(k6,k5) &
     &                         -zro(k5,k2)*zro(k6,k3)         &
     &                         +zro(k6,k2)*zro(k5,k3))        &  
     &             -zro(k5,k1)*(zkapa01(k2,k3)*zkapa10(k6,k4) &
     &                         -zro(k4,k2)*zro(k6,k3)         &
     &                         +zro(k6,k2)*zro(k4,k3))        &
     &             +zro(k6,k1)*(zkapa01(k2,k3)*zkapa10(k5,k4) &
     &                         -zro(k4,k2)*zro(k5,k3)         &
     &                         +zro(k5,k2)*zro(k4,k3))         

!        if(abs(ZRho3B).gt.0.00001) write(*,*) id123,id456,ZRho3B
!    ............................ 
!    minus sign is from the 
!    c_6c_5c_4 = -c_4c_5c_6
!    ............................ 

           ZRO3BJJ(id123,id456) = ZRO3BJJ(id123,id456) &
     &                          - ZRho3B*cg1*cg2*cg3*cg4*Znorm2/(J123+1.d0)
   
!     &           - ZRho3B*cg1*cg2*cg3*cg4*znorm/((J123+1.d0)*NFOM**2)       ! J123 is doubled  
        enddo
        enddo
        enddo
        enddo
        enddo

!     ..................................
      if(iUse_Sym.eq.1) then
        ZRO3BJJ(id456,id123) = ZRO3BJJ(id123,id456)
      endif
!     ..................................
         if(id123.eq.id456) &
     &   Trace = Trace + (J123+1.d0)*ZRO3BJJ(id123,id456)
        enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
        write(*,'(a5)',advance='no') 'done'
        time1 = omp_get_wtime()
        write(*,'(a20,f10.5)',advance='no'), '...Time elapsed:',time1-time0
        write(*,'(a14,2f15.8)') ' Trace[Rho3B]=',Trace !/Znorm2
!        write(*,'(a,i3,a,2f15.8)') &
!     &  '  Euler Angle=',iangle,'  Trace[Rho3B]=',Trace
        end


        subroutine ZRho3B_Integral(iangle,ZRO3BJJ,J,zff,ZRho3BJJ)
        USE VAPHFB_PAR
        implicit none
        integer   iangle,J
        integer   id123,id456 
        complex*16 zff
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRO3BJJ,ZRho3BJJ  

        do id123=1,Rho3B%idx123_vmax
        do id456=1,Rho3B%idx123_vmax
           ZRho3BJJ(id123,id456) = ZRho3BJJ(id123,id456)              &
     &                       +weightJ(J)*zff*ZRO3BJJ(id123,id456)  

         enddo
         enddo
       end

        subroutine Write_ME3B(ie,znorm,ZRho3BJJ)
        USE VAPHFB_PAR
        implicit none
        integer   ie
        complex*16 znorm,Trace(0:2),ZRho3B
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3BJJ  

        INTEGER id123,t123,J123,id456
!     ...............................
        Trace(0:2) = zzero
        do id123=1,Rho3B%idx123_vmax
           t123 = Rho3B%t12(id123) + Rho3B%t3(id123)
           J123 = Rho3B%J123(id123) 
        do id456=1,Rho3B%idx123_vmax

           ZRho3B = zzero

           if(abs(znorm).gt.CHOP) ZRho3B = ZRho3BJJ(id123,id456)/znorm
!     ...................................... normalization
!           ZRho3BJJ(id123,id456) = ZRho3B

           if(t123.eq.0 .and. id123.eq.id456) trace(0) = trace(0) + (J123+1.d0)*ZRho3B
           if(t123.eq.3 .and. id123.eq.id456) then
            trace(1) = trace(1) + (J123+1.d0)*ZRho3B
           endif

!           if(id123.eq.id456 .and. Rho3B%t12(id123).eq.1) trace(2) = trace(2) + (J123+1.d0)*ZRho3B 
!           pnt is not included
!           if(abs(ZRho3BJJ(id123,id456) - ZRho3BJJ(id456,id123)).gt.1.d-6) &
!     &     write(*,'(2i6,4f15.10)') id123,id456,ZRho3BJJ(id123,id456),ZRho3BJJ(id456,id123)

           if(abs(ZRho3B).gt.CHOP ) &
     &     write(ie,'(2i8,f20.10)') id123,id456, dreal(ZRho3B)
        enddo
        enddo
       write(*,*) 'Trace[Rho3B](nnn)=',Trace(0)
       write(*,*) 'Trace[Rho3B](ppp)=',Trace(1)
!       write(*,*) 'Trace[Rho3B](npt)=',Trace(2)
       end

       subroutine Write_L3B(ie,ZLambda3B)
        USE VAPHFB_PAR
        implicit none
        integer   ie
        complex*16 Trace(0:2),ZRho3B
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZLambda3B  

        INTEGER id123,t123,J123,id456
!     ...............................
        Trace(0:2) = zzero
        do id123=1,Rho3B%idx123_vmax
           t123 = Rho3B%t12(id123) + Rho3B%t3(id123)
           J123 = Rho3B%J123(id123)
        do id456=1,Rho3B%idx123_vmax

           ZRho3B = ZLambda3B(id123,id456)
!     ...................................... normalization
!           ZRho3BJJ(id123,id456) = ZRho3B

           if(t123.eq.0 .and. id123.eq.id456) trace(0) = trace(0) + (J123+1.d0)*ZRho3B
           if(t123.eq.3 .and. id123.eq.id456) then
            trace(1) = trace(1) + (J123+1.d0)*ZRho3B
           endif

           if(id123.eq.id456 .and. Rho3B%t12(id123).eq.1) trace(2) = trace(2) + (J123+1.d0)*ZRho3B 


           if(abs(ZLambda3B(id123,id456) - ZLambda3B(id456,id123)).gt.1.d-6) &
     &     write(*,'(2i6,2f20.10)') id123,id456, &
     &     dreal(ZLambda3B(id123,id456)),dreal(ZLambda3B(id456,id123))

!           if(abs(ZRho3B).gt.CHOP ) &
           if(abs(ZRho3B).gt.1.d-5 ) &
     &     write(ie,'(2i8,f20.10)') id123,id456, dreal(ZRho3B)
        enddo
        enddo
       write(*,*) 'Trace[L3B](nnn)=',Trace(0)
       write(*,*) 'Trace[L3B](ppp)=',Trace(1)
       write(*,*) 'Trace[L3B](npt)=',Trace(2)
       end



        subroutine Read_ME3B(ie,ZRho3BJJ)
        USE VAPHFB_PAR
        implicit none
        INTEGER ie,id123,id456
        complex*16 me,ZRho3B
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3BJJ  
!     ...............................
!       open(93,file='L3B_ljwang.dat',status='old')
!       open(ie,file='fort.91',status='old')


10     read(ie,*,end=20) id123,id456,me
          ZRho3BJJ(id123,id456)=me
       goto 10
20     continue
       call Check_Trace_ZRho3B(ZRho3BJJ)

!30     read(93,*,end=40) id123,id456,me
!          if(abs(ZRho3BJJ(id123,id456)-me).gt.1.d-4) &
!     &    write(500,'(2i6,4f12.8)') id123,id456,me,ZRho3BJJ(id123,id456)
!       goto 30
!40     continue
       end


        subroutine Check_Trace_ZRho3B(ZRho3BJJ)
        USE VAPHFB_PAR
        implicit none
        Complex*16 Trace(0:2),ZRho3B
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3BJJ  
        INTEGER id123,J123,t123,id456,nn,np,na
!     ...............................
        Trace(0:2) = zzero
        do id123=1,Rho3B%idx123_vmax
           J123 = Rho3B%J123(id123)
           t123 = Rho3B%t12(id123) + Rho3B%t3(id123)
!     ...................................... normalization
           ZRho3B = ZRho3BJJ(id123,id123) 
           if(t123 .eq. 0 ) trace(0) = trace(0) + (J123+1.d0)*ZRho3B
           if(t123 .eq. 3 ) trace(1) = trace(1) + (J123+1.d0)*ZRho3B
           trace(2) = trace(2) + (J123+1.d0)*ZRho3B
        enddo


        nn = Nucl%nucleon(0) - Nucl%ncore(0) 
        np = Nucl%nucleon(1) - Nucl%ncore(1)
        na = nn + np 
        if(abs(Trace(0)-nn*(nn-1)*(nn-2)).lt.1.d-8) then
           write(*,*) 'Trace[Rho3B](nnn) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](nnn)=',Trace(0)
        endif
        if(abs(Trace(1)-np*(np-1)*(np-2)).lt.1.d-8) then
           write(*,*) 'Trace[Rho3B](ppp) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](ppp)=',Trace(1)
        endif

        if(abs(Trace(2)-na*(na-1)*(na-2)).lt.1.d-8) then
           write(*,*) 'Trace[Rho3B](ttt) is OK !'
        else
           write(*,*) 'WARNING: Trace[Rho3B](ttt)=',Trace(2)
        endif
       end
!    .................................................
        subroutine ME3B_Transform(ZRho3B1,ZRho3B2)
!    .................................................
!                    --->               <----  
!    From <(12)3| | (45)6> to <(12)3| | 4(56)> 
!    .................................................
        USE VAPHFB_PAR
        implicit none
        integer   ie
        complex*16 zff,znorm,Trace(0:1),ZRho3B
        Complex*16, DIMENSION(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3B1,ZRho3B2 
        real*8  Wigner_6j
        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER k1,t1,n1,lj1,l1,j1,twom1
        INTEGER k2,t2,n2,lj2,l2,j2,twom2
        INTEGER k3,t3,n3,lj3,l3,j3,twom3
        INTEGER k4,t4,n4,lj4,l4,j4,twom4
        INTEGER k5,t5,n5,lj5,l5,j5,twom5
        INTEGER k6,t6,n6,lj6,l6,j6,twom6
        INTEGER J123,P123,P456,J12,J45,P12,P45,t12,t45
        INTEGER id123,id456,b12,a12,b45,a45,M12,M45,M123
        INTEGER tnlj1,tnlj2,tnlj4,tnlj5,tnlj6
        INTEGER SPB_l_lj,SPB_twoj_lj
        INTEGER a4,a65,b65,J65,P65,t65,id654
        INTEGER iphase_65

        ZRho3B2(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax)  = zzero
        do id123=1,Rho3B%idx123_vmax

           J123 = Rho3B%J123(id123)            ! doubled

           P12 = Rho3B%P12(id123)
           t12 = Rho3B%t12(id123)
           n3  = Rho3B%n3(id123)
           t3  = Rho3B%t3(id123)
           lj3 = Rho3B%lj3(id123)
           l3  = SPB_l_lj(lj3)*2     ! doubled l3
           j3  = SPB_twoj_lj(lj3)
           P123 = mod(l3/2+P12,2)

        do id456=1,Rho3B%idx123_vmax
           if(J123.ne.Rho3B%J123(id456)) cycle

           J45 = Rho3B%J12(id456)
           P45 = Rho3B%P12(id456)
           t45 = Rho3B%t12(id456)
           n6  = Rho3B%n3(id456)
           t6  = Rho3B%t3(id456)
           lj6 = Rho3B%lj3(id456)
! ...............

! Here l6 is zero
           l6  = SPB_l_lj(lj6)*2      ! l6 is zero ! doubled l6
           j6  = SPB_twoj_lj(lj6)
           if(lj6 .ne. (L6+j6-1)/2) stop 'Error in QNs'

           tnlj6 = Rho3B%VSPB%tnlj(lj6,n6,t6)
           P456 = mod(l6/2+P45,2)

           if(t12+t3.ne.t45+t6) cycle
           if(P123.ne.P456)     cycle

           b45   = Rho3B%VTPB%block(J45,P45,t45)
           a45   = Rho3B%a12(id456)
!      ............. quantum numbers for 1 and 2
           n4  = Rho3B%VTPB%n1(b45,a45)
           t4  = Rho3B%VTPB%t1(b45,a45)
           L4  = Rho3B%VTPB%twol1(b45,a45)  ! doubled 
           j4  = Rho3B%VTPB%j1(b45,a45)     ! doubled
           lj4 = (L4+j4-1)/2
           a4  = Rho3B%VSPB%tnlj(lj4,n4,t4)
           
           n5  = Rho3B%VTPB%n2(b45,a45)
           t5  = Rho3B%VTPB%t2(b45,a45)
           L5  = Rho3B%VTPB%twol2(b45,a45)  ! doubled 
           j5  = Rho3B%VTPB%j2(b45,a45)  ! doubled
           lj5 = (L5+j5-1)/2
           tnlj5 = Rho3B%VSPB%tnlj(lj5,n5,t5)
           
           P65 = mod((l6+l5)/2,2)
           t65 = t5+t6 

           do J65 = abs(j5-j6)/2, (j5+j6)/2  ! Not doubled
            
              b65 =  Rho3B%VTPB%block(J65,P65,t65)   
              if(t6.eq.1 .and. t5.eq.t4 .and.t5.eq.0) then ! 6(54) = p(nn)
!            ............ exchange 6 and 5 to replace pn with np
                  tnlj5 = Rho3B%VSPB%tnlj(lj6,n6,t6)
                  tnlj6 = Rho3B%VSPB%tnlj(lj5,n5,t5)
                  a65 =  Rho3B%VTPB%a(b65,tnlj6,tnlj5)  ! np
                  iphase_65 = iv(iabs(J65-(j6+j5)/2+1))
 
               else
                 a65 =  Rho3B%VTPB%a(b65,tnlj6,tnlj5) 
                 iphase_65 = 1
              endif


              id654 = Rho3B%idx(b65,a65,a4,J123)


!         if(id456.eq.261)    then
!           write(*,'(4i5)') lj6,n6,t6,Rho3B%VSPB%tnlj(lj6,n6,t6)
!           write(*,'(a,i3,a,3i3)') 'b65=',b65,'  a65=',a65,tnlj6,tnlj5 !a4,J123,id654
!               write(*,'(a,i3,a,i3,a,i3,a,4i3)') &
!     &               'J65=',J65,'t6=',t6,'t5=',t5,'t4=',t4,lj6,lj5,lj4
!         endif
 
              if(id654 .eq. -1) cycle

              ZRho3B2(id123,id456)  = ZRho3B2(id123,id456)             &
     &                               -ZRho3B1(id123,id654)             &
     &                               *dsqrt((2*J65+1.d0)*(2*J45+1.d0)) &
     &                               * iphase_65                       &
     &                             *Wigner_6j(j4,j5,2*J45,j6,J123,2*J65)


           enddo
!  ..................
         enddo !
         enddo !

        end

!    .................................................
!        subroutine ZRho3B2Lambda3B(ie,ZRho1BJ0,ZRho2BJJ,ZRho3BJJ,ZLambda3B)
        subroutine ZRho3B2Lambda3B(ie,ZRho1BJ0,ZRho3BJJ,ZLambda3B)
!    .................................................
!    ZRho1BJ0 is defined in the whole space
!    ZRhO2BJJ is defined in the pf.val space
!    ZRho3BJJ is defined in the me3b.val space 
!    .................................................
!    ATTENTION: index for quantum numbers in diff. densities might be diff.
!    .................................................
        USE VAPHFB_PAR
        implicit none
        integer    ie
        complex*16 znorm,Trace(0:1)
        Complex*16, Dimension(1:Rho3B%idx123_vmax,1:Rho3B%idx123_vmax) :: ZRho3BJJ,ZLambda3B 
!        Complex*16  ZRho2BJJ(-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
        Complex*16  ZRho1BJ0(0:1,0:HO%ljmax,0:HO%nMax,0:HO%nMax)  ! it, lj, n1,n2
        real*8  Wigner_6j
        INTEGER tt,PP,JJ,MM,bb,bm
        INTEGER k1,t1,n1,lj1,l1,j1,twom1
        INTEGER k2,t2,n2,lj2,l2,j2,twom2
        INTEGER k3,t3,n3,lj3,l3,j3,twom3
        INTEGER k4,t4,n4,lj4,l4,j4,twom4
        INTEGER k5,t5,n5,lj5,l5,j5,twom5
        INTEGER k6,t6,n6,lj6,l6,j6,twom6
        INTEGER J123,P123,P456,J12,J45,P12,P45,t12,t45
        INTEGER id123,id456,b12,a12,b45,a45,M12,M45,M123
        INTEGER tnlj1,tnlj2,tnlj4,tnlj5,tnlj6
        INTEGER SPB_l_lj,SPB_twoj_lj
        INTEGER a4

        INTEGER iphase23,J23
        INTEGER fidx1,fidx2,fidx3,fidx4,fidx5,fidx6,fb23,fa23,fa45
        INTEGER iphase64,fa64,b64,P64,t64,J64,P23,t23
        INTEGER id1,id2,id3,id4,id5,id6
        INTEGER P56,t56,iphase56,J56,fa56
        INTEGER J31,fa31,fb45,iphase31,t31,P31,fb31
        INTEGER fb12,fa12,iphase12,iphase45
!       ........ start


        ZLambda3B = ZRho3BJJ  
!      ............................ loop over id123 for lambda3B
        do id123=1,Rho3B%idx123_vmax

!          if(mod(id123,icheck) .eq. 0) write(*,*) '........ complete',id123/icheck,'%.......'
           J123 = Rho3B%J123(id123)            ! doubled

           J12 = Rho3B%J12(id123)
           P12 = Rho3B%P12(id123)
           t12 = Rho3B%t12(id123)
           t3  = Rho3B%t3(id123)
           lj3 = Rho3B%lj3(id123)
           n3  = Rho3B%n3(id123)
           l3  = SPB_l_lj(lj3)*2     ! doubled l3
           j3  = SPB_twoj_lj(lj3)    ! doubled
           P123 = mod(l3/2+P12,2)

           b12   = Rho3B%VTPB%block(J12,P12,t12)
           a12   = Rho3B%a12(id123)
!      ............. quantum numbers for 1 and 2
           n1  = Rho3B%VTPB%n1(b12,a12)
           t1  = Rho3B%VTPB%t1(b12,a12)
           L1  = Rho3B%VTPB%twol1(b12,a12)  ! doubled 
           j1  = Rho3B%VTPB%j1(b12,a12)  ! doubled
           lj1 = (L1+j1-1)/2
           n2  = Rho3B%VTPB%n2(b12,a12)
           t2  = Rho3B%VTPB%t2(b12,a12)
           L2  = Rho3B%VTPB%twol2(b12,a12)  ! doubled 
           j2  = Rho3B%VTPB%j2(b12,a12)  ! doubled
           lj2 = (L2+j2-1)/2
!      ............. check
           tnlj1 = Rho3B%VSPB%tnlj(lj1,n1,t1)
           tnlj2 = Rho3B%VSPB%tnlj(lj2,n2,t2)
!           if(a12 .ne. Rho3B%VTPB%a(b12,tnlj1,tnlj2)) &
!     &     stop 'Error in QNs for a(1,2)'
!      ............................ loop over id456 for lambda3B

        do id456=1,Rho3B%idx123_vmax
           if(J123.ne.Rho3B%J123(id456)) cycle

           J45 = Rho3B%J12(id456)
           P45 = Rho3B%P12(id456)
           t45 = Rho3B%t12(id456)
           t6  = Rho3B%t3(id456)
           lj6 = Rho3B%lj3(id456)
           n6  = Rho3B%n3(id456)
           l6  = SPB_l_lj(lj6)*2       ! doubled l6
           j6  = SPB_twoj_lj(lj6)
           P456 = mod(l6/2+P45,2)

           if(t12+t3.ne.t45+t6) cycle
           if(P123.ne.P456)     cycle

           b45   = Rho3B%VTPB%block(J45,P45,t45)
           a45   = Rho3B%a12(id456)
!      ............. quantum numbers for 1 and 2
           n4  = Rho3B%VTPB%n1(b45,a45)     ! n=0,1,2,3,...
           t4  = Rho3B%VTPB%t1(b45,a45)
           L4  = Rho3B%VTPB%twol1(b45,a45)  ! doubled 
           j4  = Rho3B%VTPB%j1(b45,a45)     ! doubled
           lj4 = (L4+j4-1)/2
           n5  = Rho3B%VTPB%n2(b45,a45)
           t5  = Rho3B%VTPB%t2(b45,a45)
           L5  = Rho3B%VTPB%twol2(b45,a45)  ! doubled 
           j5  = Rho3B%VTPB%j2(b45,a45)  ! doubled
           lj5 = (L5+j5-1)/2
           tnlj4 = Rho3B%VSPB%tnlj(lj4,n4,t4)
           tnlj5 = Rho3B%VSPB%tnlj(lj5,n5,t5)

!       here

        t23 = t2 + t3
        P23 = mod((l2+l3)/2,2) 

        t64 = t6+t4
        P64 = mod((l6+l4)/2,2)

        t31 = t1+t3
        P31 = mod((l1+l3)/2,2)

        t56 = t5+t6
        P56 = mod((l5+l6)/2,2)

!        if(id123.ge.66 .and. id123 .le.76 .and. id456.eq.id123) then
!        if(id123.eq.1896 .and. id456.eq.554) then         
!         write(*,'(a,i5,a,i5)') 'id123=',id123,' id456=',id456
!         write(*,'(a,3i5)') '1',t1,n1,lj1
!         write(*,'(a,3i5)') '2',t2,n2,lj2
!         write(*,'(a,3i5)') '3',t3,n3,lj3
!         write(*,'(a,3i5)') '4',t4,n4,lj4
!         write(*,'(a,3i5)') '5',t5,n5,lj5
!         write(*,'(a,3i5)') '6',t6,n6,lj6
!         write(*,'(a,i5,a,i5,a,i5)') 'J12=',J12,' J45=',J45,' J123=',J123
!         stop
!        endif
       
!  ..............................................
!  ..............................................
!  ..............................................
!  ..............................................
!  ..............................................
!
!  .... T1 
        if(lj1.eq.lj6 .and. t1.eq.t6 .and. P45.eq.P23 .and. t23.eq.t45) then      

          J23 = J45
          if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else  ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb23 = TPB%block(J23,P23,t23) 
           fa23 = TPB%a(fb23,fidx2,fidx3) 
           fa45 = TPB%a(fb23,fidx4,fidx5) 

!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J45+(j2+j3)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *iphase23*iphase45*Dens%ZRho2BJJ(fa23,fa45,fb23)      ! n=1,2,..
          endif

!  .... T2 

        if(lj1.eq.lj5 .and. t1.eq.t5 .and. P64.eq.P23 .and. t64.eq.t23) then

         do J23 = abs(j2-j3)/2, (j2+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(23)(64), if 2=3, or 4=6, J is only even
            if(t2.eq.t3 .and. lj2.eq.lj3 .and. iv(J23).eq.-1) cycle
            if(t4.eq.t6 .and. lj4.eq.lj6 .and. iv(J23).eq.-1) cycle
            if(J23.lt.abs(j6-j4)/2 .or. J23.gt.(j6+j4)/2) cycle

           J64 = J23
           if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else    ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t23.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif 
!      ................(2,3)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t23.ne.1.and.id6.gt.id4) then
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j6+j4)/2+1)
            endif
 

           fb23 = TPB%block(J23,P23,t23)
           fa23 = TPB%a(fb23,fidx2,fidx3)
           fa64 = TPB%a(fb23,fidx6,fidx4)

           if(fa23.eq.-1 .or. fa64.eq.-1) then
!           write(*,'(10i6)') J23,P23,t23,fb23,fidx2,fidx3,fa23,fidx6,fidx4,fa64 
!           write(*,*) fb23,fa23,fa64 
!           write(*,*) 't',t2,t3,t6,t4 
!           write(*,*) 'n',n2,n3,n6,n4 
!           write(*,*) 'l',l2,l3,l6,l4 
           write(*,*) 'j',j2,j3,j6,j4 
           stop
           endif  
           
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv(J45+J23+(j1+j2+j3+j4)/2) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J23+1) &
     &                 *Wigner_6j(j4,j1,2*J45,J123,j6,2*J23) &
     &                 *Wigner_6j(j2,j3,2*J23,J123,j1,2*J12) &
     &                 *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                 *iphase23*iphase64*Dens%ZRho2BJJ(fa23,fa64,fb23)      ! n=1,2,..
          enddo
        endif

!  .... T3 

        if(lj1.eq.lj4 .and. t1.eq.t4 .and. &
     &     mod((l6+l5)/2,2) .eq. P23 .and. t23.eq.t5+t6) then

         do J23 = abs(j2-j3)/2, (j2+j3)/2 ! NOT Doubled

!        .......................................... Rho2B(23)(64), if 2=3, or 4=6, J is only even
            if(t2.eq.t3 .and. lj2.eq.lj3 .and. iv(J23).eq.-1) cycle
            if(t5.eq.t6 .and. lj5.eq.lj6 .and. iv(J23).eq.-1) cycle
            if(J23.lt.abs(j6-j5)/2 .or. J23.gt.(j6+j5)/2) cycle

           J56 = J23
           if(t2.eq.1 .and.t3.eq.0) then
!      .................................. (23) = (pn)
           fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = iv(J23-(j2+j3)/2+1)
          else    ! nn,np,pp
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase23 = 1
          endif

           if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1
          endif

            id2 = n2+1 + (HO%NMax)*lj2
            id3 = n3+1 + (HO%NMax)*lj3
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t23.ne.1.and.id2.gt.id3) then
              fidx2 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase23 = iv(J23-(j2+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t23.ne.1.and.id5.gt.id6) then
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j6+j5)/2+1)
            endif


           fb23 = TPB%block(J23,P23,t23)
           fa23 = TPB%a(fb23,fidx2,fidx3)
           fa56 = TPB%a(fb23,fidx5,fidx6)

           if(fa23.eq.-1 .or. fa56.eq.-1) then
           write(*,*) 'j',j2,j3,j5,j6
           stop
           endif

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j2+j3+j5+j6)/2) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J23+1) &
     &                 *Wigner_6j(j5,j1,2*J45,J123,j6,2*J23) &
     &                 *Wigner_6j(j2,j3,2*J23,J123,j1,2*J12) &
     &                 *ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                 *iphase23*iphase56*Dens%ZRho2BJJ(fa23,fa56,fb23)      ! n=1,2,..
          enddo
        endif

4     continue

!  .... T4 

        if(lj2.eq.lj6 .and. t2.eq.t6 .and. P45 .eq. P31 .and. t31.eq.t45) then

          J31 = J45
          if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
          else  ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n3 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif

           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id3.gt.id1) then
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb45 = TPB%block(J45,P45,t45)
           fa31 = TPB%a(fb45,fidx3,fidx1)
           fa45 = TPB%a(fb45,fidx4,fidx5)

!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv((j1+j2)/2-J12+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *iphase31*iphase45*Dens%ZRho2BJJ(fa31,fa45,fb45)      ! n=1,2,..
          endif



!  .... T5 

        if(lj2.eq.lj5 .and. t2.eq.t5 .and. P64.eq.P31 .and. t31.eq.t64) then

         do J31 = abs(j1-j3)/2, (j1+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(31)(64), if 2=3, or 4=6, J is only even
            if(t1.eq.t3 .and. lj1.eq.lj3 .and. iv(J31).eq.-1) cycle
            if(t4.eq.t6 .and. lj4.eq.lj6 .and. iv(J31).eq.-1) cycle
            if(J31.lt.abs(j6-j4)/2 .or. J31.gt.(j6+j4)/2) cycle

           J64 = J31
           if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
          else    ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif

           if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t31.ne.1.and.id3.gt.id1) then
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t31.ne.1.and.id6.gt.id4) then
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j6+j4)/2+1)
            endif


           fb31 = TPB%block(J31,P31,t31)
           fa31 = TPB%a(fb31,fidx3,fidx1)
           fa64 = TPB%a(fb31,fidx6,fidx4)

           if(fa31.eq.-1 .or. fa64.eq.-1) then
           write(*,*) 'j',j3,j1,j6,j4
           stop
           endif

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j4+j1)/2+J12+J45+1) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J31+1) &
     &                 *Wigner_6j(j4,j2,2*J45,J123,j6,2*J31) &
     &                 *Wigner_6j(j1,j3,2*J31,J123,j2,2*J12) &
     &                 *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                 *iphase31*iphase64*Dens%ZRho2BJJ(fa31,fa64,fb31)      ! n=1,2,..
          enddo
        endif


!  .... T6 
        if(lj2.eq.lj4 .and. t2.eq.t4 .and. P56 .eq. P31 .and. t31.eq.t56) then

         do J31 = abs(j1-j3)/2, (j1+j3)/2 ! NOT Doubled
!        .......................................... Rho2B(31)(64), if 2=3, or 4=6, J is only even
            if(t1.eq.t3 .and. lj1.eq.lj3 .and. iv(J31).eq.-1) cycle
            if(t5.eq.t6 .and. lj5.eq.lj6 .and. iv(J31).eq.-1) cycle
            if(J31.lt.abs(j6-j5)/2 .or. J31.gt.(j6+j5)/2) cycle

           J56 = J31
           if(t3.eq.1 .and.t1.eq.0) then
!      .................................. (31) = (pn)
           fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = iv(J31-(j1+j3)/2+1)
          else    ! nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           fidx3 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase31 = 1
          endif

           if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1
          endif

            id3 = n3+1 + (HO%NMax)*lj3
            id1 = n1+1 + (HO%NMax)*lj1
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(3,1)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t31.ne.1.and.id3.gt.id1) then
              fidx3 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj3,n3+1,t3)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase31 = iv(J31-(j1+j3)/2+1)
            endif
!      ................(5,6)=nn, or pp, in which case, only id6 =< id4 is stored in TPB
            if(t31.ne.1.and.id5.gt.id6) then
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j6+j5)/2+1)
            endif


           fb31 = TPB%block(J31,P31,t31)
           fa31 = TPB%a(fb31,fidx3,fidx1)
           fa56 = TPB%a(fb31,fidx5,fidx6)

           if(fa31.eq.-1 .or. fa56.eq.-1) then
           write(*,*) 'j',j3,j1,j5,j6
           stop
           endif

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)   &
     &                         - iv((j1+j2+j5+j6)/2+J12+J31) &
     &                 *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)*(2*J31+1) &
     &                 *Wigner_6j(j5,j2,2*J45,J123,j6,2*J31) &
     &                 *Wigner_6j(j1,j3,2*J31,J123,j2,2*J12) &
     &                 *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                 *iphase31*iphase56*Dens%ZRho2BJJ(fa31,fa56,fb31)      ! n=1,2,..
          enddo
        endif



!  .... T7 

        if(lj3.eq.lj6 .and. t3.eq.t6 .and. J12.eq.J45 .and. P12.eq.P45 .and. t12.eq.t45) then


!         .............(1,2) only contains nn,np,pp
           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

!         .............(4,5) only contains nn,np,pp
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n5 starts from 0, while n in FSPB starts from 1,.. 
           iphase45 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id4 = n4+1 + (HO%NMax)*lj4
            id5 = n5+1 + (HO%NMax)*lj5
!      ................(1,2)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id1.gt.id2) then
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j1+j2)/2+1)
            endif
!      ................(4,5)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t45.ne.1.and.id4.gt.id5) then
              fidx4 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx5 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase45 = iv(J45-(j4+j5)/2+1)
            endif

           fb45 = TPB%block(J45,P45,t45)
           fa12 = TPB%a(fb45,fidx1,fidx2)
           fa45 = TPB%a(fb45,fidx4,fidx5)

!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - ZRho1BJ0(t3,lj3,n3,n6)               & ! n=0,1,2,..
     &                          *iphase12*iphase45*Dens%ZRho2BJJ(fa12,fa45,fb45)      ! n=1,2,..
          endif


8      continue        
!  .... T8 
        if(lj3.eq.lj5 .and. t3.eq.t5 .and. P12.eq.P64  .and. t12.eq.t64) then

          J64 = J12
          if(t6.eq.1 .and.t4.eq.0) then
!      .................................. (64) = (pn)
           fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = iv(J64-(j6+j4)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx4 = FSPB%tnlj(lj4,n4+1,t4)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase64 = 1
          endif

           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id6 = n6+1 + (HO%NMax)*lj6
            id4 = n4+1 + (HO%NMax)*lj4
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id1.gt.id2) then
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j2+j1)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id6.gt.id4) then
              fidx4 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx6 = FSPB%tnlj(lj4,n4+1,t4)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase64 = iv(J64-(j4+j6)/2+1)
            endif

!           do J23=abs(j2-j3)/2, (j2+j3)/2   ! NOT doubled 
           fb12 = TPB%block(J12,P12,t12)
           fa12 = TPB%a(fb12,fidx1,fidx2)
           fa64 = TPB%a(fb12,fidx6,fidx4)

!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J45+(j3+j4)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j4,j3,2*J45,J123,j6,2*J12) &
     &                          *ZRho1BJ0(t3,lj3,n3,n5)               & ! n=0,1,2,..
     &                          *iphase12*iphase64*Dens%ZRho2BJJ(fa12,fa64,fb12)      ! n=1,2,..
          endif


9      continue        


!  .... T9 
        if(lj3.eq.lj4 .and. t3.eq.t4 .and. P12.eq.P56  .and. t12.eq.t56) then

          J56 = J12
          if(t5.eq.1 .and.t6.eq.0) then
!      .................................. (56) = (pn)
           fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = iv(J56-(j6+j5)/2+1)
          else  ! nn,np,pp
           fidx6 = FSPB%tnlj(lj6,n6+1,t6)  ! n6 starts from 0, while n in FSPB starts from 1,.. 
           fidx5 = FSPB%tnlj(lj5,n5+1,t5)  ! n4 starts from 0, while n in FSPB starts from 1,.. 
           iphase56 = 1
          endif

           fidx1 = FSPB%tnlj(lj1,n1+1,t1)  ! n1 starts from 0, while n in FSPB starts from 1,.. 
           fidx2 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
           iphase12 = 1

            id1 = n1+1 + (HO%NMax)*lj1
            id2 = n2+1 + (HO%NMax)*lj2
            id5 = n5+1 + (HO%NMax)*lj5
            id6 = n6+1 + (HO%NMax)*lj6
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id1.gt.id2) then
              fidx2 = FSPB%tnlj(lj1,n1+1,t1)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx1 = FSPB%tnlj(lj2,n2+1,t2)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase12 = iv(J12-(j2+j1)/2+1)
            endif
!      ................(2,3)=nn, or pp, in which case, only id2 =< id3 is stored in TOB
            if(t12.ne.1.and.id5.gt.id6) then
              fidx5 = FSPB%tnlj(lj6,n6+1,t6)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              fidx6 = FSPB%tnlj(lj5,n5+1,t5)  ! n2 starts from 0, while n in FSPB starts from 1,.. 
              iphase56 = iv(J56-(j5+j6)/2+1)
            endif

!           do J23=abs(j2-j3)/2, (j2+j3)/2   ! NOT doubled 
           fb12 = TPB%block(J12,P12,t12)
           fa12 = TPB%a(fb12,fidx1,fidx2)
           fa56 = TPB%a(fb12,fidx5,fidx6)

!           write(*,'(10i6)') J23,P45,t45,fb23,fidx2,fidx3,fa23,fidx4,fidx5,fa45 
!           write(*,*) fb23,fa23,fa45 
           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         - iv(J12+(j5+j6)/2+1)*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j5,j3,2*J45,J123,j6,2*J12) &
     &                          *ZRho1BJ0(t3,lj3,n3,n4)               & ! n=0,1,2,..
     &                          *iphase12*iphase56*Dens%ZRho2BJJ(fa12,fa56,fb12)      ! n=1,2,..
          endif



10      continue        
!  .... T10 
        if( lj1.eq.lj6 .and. t1.eq.t6 .and. lj2.eq.lj5 .and. t2.eq.t5 .and. lj3.eq.lj4 .and. t3.eq.t4) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456) &
     &                         +2*dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0) &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n4)                 ! n=0,1,2,..
          endif


!  .... T11 
        if( lj1.eq.lj6 .and. t1.eq.t6 .and. lj2.eq.lj4 .and. t2.eq.t4 .and. lj3.eq.lj5 .and. t3.eq.t5) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         +2*iv(J45+(j2+j3)/2+1)                 &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j1,j2,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n5)                 ! n=0,1,2,..
          endif


!  .... T12 
        if( lj1.eq.lj5 .and. t1.eq.t5 .and. lj2.eq.lj6 .and. t2.eq.t6 .and. lj3.eq.lj4 .and. t3.eq.t4) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         +2*iv((j1+j2)/2-J12+1)                 &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n4)                 ! n=0,1,2,..
          endif

!  .... T13 
        if( lj1.eq.lj5 .and. t1.eq.t5 .and. lj2.eq.lj4 .and. t2.eq.t4 .and. lj3.eq.lj6 .and. t3.eq.t6) then

        if(J12.eq.J45) ZLambda3B(id123,id456) = ZLambda3B(id123,id456)&
     &                         -2*iv((j1+j2)/2-J12)                   &
     &                          *ZRho1BJ0(t1,lj1,n1,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n6)                 ! n=0,1,2,..
          endif

14      continue        
!  .... T14 
        if( lj1.eq.lj4 .and. t1.eq.t4 .and. lj2.eq.lj6 .and. t2.eq.t6 .and. lj3.eq.lj5 .and. t3.eq.t5) then

           ZLambda3B(id123,id456) = ZLambda3B(id123,id456)            &
     &                         -2*iv((j3+j2)/2-J12+J45)               &
     &                          *dsqrt(2*J12+1.d0)*dsqrt(2*J45+1.d0)  &
     &                          *Wigner_6j(j2,j1,2*J12,j3,J123,2*J45) &
     &                          *ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n6)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n5)                 ! n=0,1,2,..
          endif

15      continue        
!  .... T15 
        if( lj1.eq.lj4 .and. t1.eq.t4 .and. lj2.eq.lj5 .and. t2.eq.t5 .and. lj3.eq.lj6 .and. t3.eq.t6) then

        if(J12.eq.J45) ZLambda3B(id123,id456) = ZLambda3B(id123,id456)&
     &                        +2*ZRho1BJ0(t1,lj1,n1,n4)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t2,lj2,n2,n5)               & ! n=0,1,2,..
     &                          *ZRho1BJ0(t3,lj3,n3,n6)                 ! n=0,1,2,..
          endif
!  ..............................................
           enddo
           enddo
        end
