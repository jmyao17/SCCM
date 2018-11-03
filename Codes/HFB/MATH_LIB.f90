      subroutine gfv() 

!=======================================================================
!
!     Calculates sign, sqrt, factorials, etc. of integers and half int.
!
!     iv(n)  =  (-1)**n
!     sq(n)  =  sqrt(n)
!     sqi(n) =  1/sqrt(n)
!     sqh(n) =  sqrt(n+1/2)
!     shi(n) =  1/sqrt(n+1/2)
!     fak(n) =  n!
!     ibc(m,n) = m!/(n!(m-n)!)  
!     fad(n) =  (2*n+1)!!
!     fdi(n) =  1/(2*n+1)!!
!     fi(n)  =  1/n!
!     wf(n)  =  sqrt(n!)
!     wfi(n) =  1/sqrt(n!)
!     wfd(n) =  sqrt((2*n+1)!!)
!     gm2(n) =  gamma(n+1/2)
!     gmi(n) =  1/gamma(n+1/2)
!     wg(n)  =  sqrt(gamma(n+1/2))
!     wgi(n) =  1/sqrt(gamma(n+1/2))
!
!-----------------------------------------------------------------------
      USE VAPHFB_PAR 
      Implicit none
      INTEGER i,j,m,n
!      INTEGER iv(-igfv:igfv),ibc(0:igfvbc,0:igfvbc)
!      Real*8  sq(0:igfv),sqi(0:igfv),sqh(0:igfv),shi(0:igfv)
!      Real*8  fak(0:igfv),fad(0:igfv),fi(0:igfv),fdi(0:igfv),wf(0:igfv),wfi(0:igfv)
!      Real*8  wfd(0:igfv),gm2(0:igfv),gmi(0:igfv),wg(0:igfv),wgi(0:igfv)

      mit(1)  = 0  ! n
      mit(-1) = 1  ! p
      mip(1)  = 1 ! +
      mip(-1) = 2 ! -

      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  one
      fad(0) =  one
      fi(0)  =  one
      fdi(0) =  one
      wf(0)  =  one
      wfi(0) =  one
      wfd(0)=  one
!     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
      do i = 1,igfv
         iv(i)         = -iv(i-1)
         iv(-igfv+i-1) = -iv(i)
         sq(i)  = dsqrt(dfloat(i))
         sqi(i) = one/sq(i)
         sqh(i) = sqrt(i+half)
         shi(i) = one/sqh(i)
         fak(i) = i*fak(i-1)
         fak(-i) = fak(i)
         fad(i) = (2*i+1)*fad(i-1)
         fi(i)  = one/fak(i)
         fdi(i) = one/fad(i)
         wf(i)  = sq(i)*wf(i-1)
         wfi(i) = one/wf(i)
         wfd(i) = sqrt(fad(i))
         gm2(i) = (i-half)*gm2(i-1)
         gmi(i) = one/gm2(i)
         wg(i)  = sqh(i-1)*wg(i-1)
         wgi(i) = one/wg(i)
      enddo
      ibc(0,0)= one
      do m=1,igfv !bc
        do n=0,m
           ibc(m,n)=fak(m)/(fak(n)*fak(m-n))
        enddo
      enddo
      end


       double precision function rnla(n1,nl1,nl,n2,nl2)

!=======================================================================
!      calculate <n1,l1| (r/b)^l |n2,l2>
!      where n1=1,2,3,4,...
!======================================================================c
!     calculates the radial functions for the spherical oscillator
!     the wave function R_nl(r) of the spherical oscillator are: 
!     phi(r,Omega) = b^(-3/2) * R_nl(r) * Y_ljm(Omega) 
!     R_nl(r) = N_nl * r**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
!     N_nl    = sqrt(2 * (n-1)!/(n+l-1/2)!)     and    x=r/b
!     n=1,2,3,...    
!     R_nl is normalized in such way that the norm integral reads
!     \int dr r**2 R_nl(r)^2 = 1 
!     However, in the 
!     Ref: S.G. Nilsson, Mat-fys Medd 29 (1955)16,1-69
!     the radial quantum number n=n-1, i.e.,
!      R_nl(r) = N_nl * r**l * L^(l+1/2)_n(x*x) * exp(-x*x/2), n=0,1,2,...
!      one should pay attention to
!-----------------------------------------------------------------------

      USE VAPHFB_PAR 
      IMPLICIT NONE

      Integer it,imu1,imu2      
      Integer n1,n2,nn1,nn2,nl,nl1,nl2,is,ismin,ismax
        
      Real*8  fac1,facs

      rnla = 0.d0
      nn1  = n1-1
      nn2  = n2-1

!       if(nl1.gt.(nl2+nl).or.nl1.lt.abs(nl2-nl))      return  ! <l1|| Y_L || l2 >
!       if((2*nn1+nl1).gt.(2*nn2+nl2+nl).or.           &
!     &   (2*nn1+nl1).lt.(2*nn2+nl2-nl))               return
!       if(iv(nl1-nl2+nl).ne.1.or.iv(nl2-nl1+nl).ne.1) return  ! parity


       imu2  = (nl1-nl2+nl)/2   !nu
       imu1  = (nl2-nl1+nl)/2   !nu'
       it    = (nl2+nl1+nl)/2   !t-1/2

!      .............. gm2(n) = Gamma(n+1/2)
       fac1  = fak(nn1)*fak(nn2)/(gm2(nn2+nl2+1)*gm2(nn1+nl1+1))
       facs  = sqrt(fac1)*fak(imu1)*fak(imu2)
       ismin = max(nn1-imu1,nn2-imu2)
       ismax = min(nn1,nn2)

       do is=ismin,ismax

         if(is.lt.0) goto 10
         if((nn1-is).lt.0) goto 10
         if((nn2-is).lt.0) goto 10
         if((is+imu2-nn2).lt.0) goto 10
         if((is+imu1-nn1).lt.0) goto 10

         rnla = rnla + gm2(it+is+1)*facs                          &
     &                /(fak(is)*fak(nn2-is)*fak(nn1-is)           &
     &                *fak(is+imu2-nn2)*fak(is+imu1-nn1)) 

  10     enddo !is


       rnla=rnla*iv(nn1+nn2)

      return
!-----end rnl
      end


!=======================================================================

      subroutine sdiag(nmax,n,a,d,x,e,is)

!=======================================================================
!
!     A   matrix to be diagonalized
!     D   eigenvalues    
!     X   eigenvectors
!     E   auxiliary field
!     IS = 1  eigenvalues are ordered and major component of X is positiv
!          0  eigenvalues are not ordered            
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
!
      dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)
!
      data tol,eps/1.e-32,1.e-10/
!
      if (n.eq.1) then
         d(1)=a(1,1)
         x(1,1)=1.
         return
      endif

      do 10 i=1,n
      do 10 j=1,i
   10    x(i,j)=a(i,j)
!
!cc   householder-reduktion
      i=n
   15 if (i-2) 200,20,20
   20 l=i-2
      f=x(i,i-1)
      g=f
      h=0
      if (l) 31,31,32
   32 do 30 k=1,l
   30 h=h+x(i,k)*x(i,k)
   31 s=h+f*f
      if (s-tol) 33,34,34
   33 h=0
      goto 100
   34 if (h) 100,100,40
   40 l=l+1
      g= dsqrt(s)
      if (f.ge.0.) g=-g
      h=s-f*g
      hi=1.d0/h
      x(i,i-1)=f-g
      f=0.0
      if (l) 51,51,52
   52 do 50 j=1,l

      x(j,i)=x(i,j)*hi
      s=0.0
      do 55 k=1,j
   55 s=s+x(j,k)*x(i,k)
      j1=j+1
      if (l-j1) 57,58,58
   58 do 59 k=j1,l
   59 s=s+x(k,j)*x(i,k)
   57 e(j)=s*hi
   50 f=f+s*x(j,i)
   51 f=f*hi*.5d0
!                                    
      if (l) 100,100,62
   62 do 60 j=1,l
      s=x(i,j)
      e(j)=e(j)-f*s
      p=e(j)
      do 65 k=1,j
   65 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
   60 continue
  100 continue
      d(i)=h
      e(i-1)=g
      i=i-1
      goto 15
!cc   Bereitstellen der Transformationmatrix 
  200 d(1)=0.0
      e(n)=0.0
      b=0.0
      f=0.0
      do 210 i=1,n
      l=i-1
      if (d(i).eq.0.) goto 221
      if (l) 221,221,222
  222 do 220 j=1,l
      s=0.0
      do 225 k=1,l
  225 s=s+x(i,k)*x(k,j)
      do 226 k=1,l
  226 x(k,j)=x(k,j)-s*x(k,i)
  220 continue
  221 d(i)=x(i,i)
      x(i,i)=1
      if (l) 210,210,232
  232 do 230 j=1,l
      x(i,j)=0.0
  230 x(j,i)=0.0
  210 continue
!
!cc   Diagonalisieren der Tri-Diagonal-Matrix
      DO 300 L=1,N
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h
!cc   Test fuer Splitting        
      do 310 j=l,n
      if ( abs(e(j)).le.b) goto 320
  310 continue
!
!cc   test fuer konvergenz    
  320 if (j.eq.l) goto 300
  340 p=(d(l+1)-d(l))/(2*e(l))
      r= dsqrt(p*p+1.d0)
      pr=p+r
      if (p.lt.0.) pr=p-r
      h=d(l)-e(l)/pr
      do 350 i=l,n
  350 d(i)=d(i)-h
      f=f+h
!
!cc   QR-transformation          
      p=d(j)
      c=1.d0
      s=0.0
      i=j
  360 i=i-1
      if (i.lt.l) goto 362
      g=c*e(i)
      h=c*p
      if ( abs(p)- abs(e(i))) 363,364,364
  364 c=e(i)/p
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.d0/r
      goto 365
  363 c=p/e(i)
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r
      s=1.d0/r
      c=c/r
  365 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 368 k=1,n
         h=x(k,i+1)
         x(k,i+1)=x(k,i)*s+h*c
  368    x(k,i)=x(k,i)*c-h*s
      goto 360
  362 e(l)=s*p
      d(l)=c*p
      if ( abs(e(l)).gt.b) goto 340
!
!cc   konvergenz      
  300 d(l)=d(l)+f
!
      if (is.eq.0) return
!cc   ordnen der eigenwerte    
      do 400 i=1,n
      k=i
      p=d(i)
      j1=i+1
      if (j1-n) 401,401,400
  401 do 410 j=j1,n
      if (d(j).ge.p) goto 410
      k=j
      p=d(j)
  410 continue
  420 if (k.eq.i) goto 400
      d(k)=d(i)
      d(i)=p
      do 425 j=1,n
      p=x(j,i)
      x(j,i)=x(j,k)
  425 x(j,k)=p
  400 continue
!                 
!     signum
      do k = 1,n
         s = 0.0d0
         do i = 1,n
            h = abs(x(i,k))
            if (h.gt.s) then
               s  = h
               im = i
            endif
         enddo   ! i
         if (x(im,k).lt.0.0d0) then
            do i = 1,n
               x(i,k) = - x(i,k)
            enddo
         endif
      enddo   ! k
! 
      return
!-end-SDIAG
      end


 
        subroutine progress(j)
        implicit none
        integer(kind=4)::j,k
        character(len=17):: bar="???% |          |"
        write(unit=bar(1:3),fmt="(i3)")  j ! 10*j
        do k=1, j
           bar(6+k:6+k)="*"
        enddo
! print the progress bar.
        write(unit=6,fmt="(a1,a1,a17)") '+',char(13), bar
        return
        end subroutine progress

!======================================================================c
      subroutine degen(na,n,ea,dd,bb,eb,eps,zz,z)
!
!======================================================================c
!
!     EA    is a set of partially degenerate eigenvalues of some matrix 
!     DD    are the corresponding eigenvectors DD. 
!     BB    is a matrix, which is diagonalized in the subspaces of
!           degenerate eigenvales EA.
!     EB    contains the eigenvalues of BB in these subspaces.
!     EPS   determines, to which accuracy the eigenvalues EA are
!           considered to be degenerate
!     ZZ,Z  are auxiliary arrays
!
!----------------------------------------------------------------------c
      use vaphfb_par
      implicit real*8 (a-h,o-z) 
!
      dimension bb(na,n),dd(na,n),ea(n),eb(n)
      dimension zz(na,n),z(n)
!
!     common /mathco/ zero,one,two,half,third,pi
!      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka
! 
!---- check for degeneracies
      k1 = 0
      do i1 = 1,n
         k1 = k1 + 1
         if (k1.ge.n) goto 20
         do k2 = k1+1,n
            if (ea(k2)-ea(k1).gt.eps) goto 10
         enddo
   10    me = k2 - 1
         ma = k1
!
!----    diagonalize together with bb 
         if (me.gt.k1) then
            m0 = ma - 1
            mm = me - m0

            do m1 = ma,me
               do k = 1,n
                  s = zero
                  do i = 1,n  
                     s = s + dd(i,m1) * bb(i,k)
                  enddo   ! i
                  z(k) = s
               enddo   ! k
               do m2 = ma,me
                  s = zero
                  do k = 1,n
                     s = s + z(k) * dd(k,m2)
                  enddo   ! k
                  zz(m1-m0,m2-m0) = s
               enddo   ! m2
            enddo   ! m1
            call sdiag(na,mm,zz,eb(ma),zz,z,+1)
            do i = 1,n
               do m = ma,me
                  s = zero
                  do l = ma,me 
                     s = s + dd(i,l) * zz(l-m0,m-m0)
                  enddo   ! l
                  z(m) = s
               enddo   ! m
               do m = ma,me
                  dd(i,m) = z(m)
               enddo   ! m
            enddo   ! i
            k1 = me 
         endif
      enddo   ! i1
!
   20 return
!-end-DEGEN
      end
