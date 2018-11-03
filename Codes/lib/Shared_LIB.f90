      subroutine StripSpaces(string)
      character(len=*) :: string
      integer :: stringLen
      integer :: last, actual

      stringLen = len (string)
      last = 1
      actual = 1

      do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
      end do
      end subroutine
!=======================================================================
      subroutine nucleus(is,npro,te)
!=======================================================================
!     is = 1 determines the symbol for a given proton number npro
!          2 determines the proton number for a given symbol te
!-----------------------------------------------------------------------
      PARAMETER (MAXZ=140)
      CHARACTER TE*2,T*(2*MAXZ+2)
!
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
!
! ... Rf is called also as Ku (kurchatovium)
! ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
! ... Ha is called also as Db (Dubnium)
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong'
         te = t(2*npro+1:2*npro+2)
         return
      else
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
               if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
      stop
      END

        integer function mitf(mt)
        integer mt

        if(mt .eq. -1) mitf = 1   ! p
        if(mt .eq.  1) mitf = 0   ! n
        if(abs(mt).ne.1) stop 'mt should be +1 or -1'

        return
        end

        subroutine print_matrix(LDA,A_matrix)
        implicit none
        integer LDA,ii,jj
        real*8  A_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,f12.8)') jj,ii,A_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_zmatrix(LDA,ZA_matrix)
        implicit none
        integer     LDA,ii,jj
        complex*16  ZA_matrix(LDA,LDA)

        do ii=1,LDA
        do jj=1,LDA
           write(911,'(2i6,2e12.3)') jj,ii,ZA_matrix(jj,ii)
        enddo
        enddo
        stop
        END SUBROUTINE

        subroutine print_array(LDA,V,istop)
        implicit none
        integer     LDA,ii,istop
        real*8 V(LDA)

        do ii=1,LDA
           write(911,'(i6,e12.3)') ii,V(ii)
        enddo
        if(istop .ne.0) stop

        END SUBROUTINE

!______________________________________________________________________________
      subroutine clingd(ma,mx,n,m,a,x,d,ifl)
!======================================================================c
!     solves the system A*X=B, where B is at the beginning on X
!     it will be overwritten later on, d is the determinant
!----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
!
      COMPLEX*16 A,X,D,CP,CQ
!
      DIMENSION  A(MA,N),X(MX,M)

      DATA CC/4H----/,TOLLIM/1.E-12/
      IFL=1
      P=0.D0
      DO 10 I=1,N
      Q=0.D0
      DO 20 J=1,N
 20   Q=Q+CDABS(A(I,J))
      IF (Q.GT.P)   P=Q
 10   CONTINUE

      TOL=TOLLIM*P
      D=CMPLX(1.D0,0.D0)

      DO 30 K=1,N
      P=0.D0
      DO 40 J=K,N
      Q=CDABS(A(J,K))
      IF (Q.LT.P) GOTO 40
      P=Q
      I=J
 40   CONTINUE
      IF (P.GT.TOL) GOTO 70
      WRITE (6,200) (CC,J=1,22),TOL,I,K,A(I,K),(CC,J=1,22)
  200 FORMAT (/1X,22A4/' *****  FAIL IN LINGM, TOLERANZ =',E10.4,  &
     & ' FROM A(',I3,',',I3,') IS ',E10.4/1X,22A4)
      IFL=-1
      RETURN

   70 CP=1./A(I,K)
      IF (I.EQ.K) GOTO 90
      D=-D

      DO 81 L=1,M
      CQ=X(I,L)
      X(I,L)=X(K,L)
   81 X(K,L)=CQ
      DO 80 L=K,N
      CQ=A(I,L)
      A(I,L)=A(K,L)
   80 A(K,L)=CQ
   90 D=D*A(K,K)
      IF (K.EQ.N) GOTO 1
      K1=K+1
      DO 120 I=K1,N
      CQ=A(I,K)*CP
      DO 106 L=1,M
  106 X(I,L)=X(I,L)-CQ*X(K,L)
      DO 120 L=K1,N
  120 A(I,L)=A(I,L)-CQ*A(K,L)
   30 CONTINUE

    1 DO 126 L=1,M
  126 X(N,L)=X(N,L)*CP
      IF (N.EQ.1) RETURN
      N1=N-1
      DO 140 K=1,N1
      CP=1./A(N-K,N-K)
      DO 140 L=1,M
      CQ=X(N-K,L)
      DO 141 I=1,K
  141 CQ=CQ-A(N-K,N+1-I)*X(N+1-I,L)
  140 X(N-K,L)=CQ*CP

      return
      end


!_______________________________________________________________________
      subroutine sdiag(nmax,n,a,d,x,e,is)

!......................................................................
!
!     A   matrix to be diagonalized
!     D   eigenvalues    
!     X   eigenvectors
!     E   auxiliary field
!     IS = 1  eigenvalues are ordered and major component of X is positiv
!          0  eigenvalues are not ordered            
!......................................................................
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
!
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
      DO 300 L=1,N
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h
!
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



!=======================================================================
      subroutine sdiag2(nmax,n,a,d,x,e,is)
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
!==================================================================
      double precision function djmk(j,m,k,cosbet,is)
!================================================================== 
!     Calculates the Wigner-Functions  d_mm'^j (beta) 
!     IS = 0: integer values for j = J,  m = M, m' = K
!
!     IS = 1: half integer valus for j = J-1/2 , m = M-1/2, m' = K-1/2
!
!     COSBET = cos(beta)
!
!     The Wigner-Functions and the phases are defined as in Edmonds.
!
!-----------------------------------------------------------------------
      USE VAPHFB_Par
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      djmk = 0.d0

      if (iabs(m).gt.j.or.iabs(k).gt.j) return

      if (abs(cosbet).gt.1.d0)   &
     &   stop 'in DJMK: cos(beta) is larger than one'

      if (cosbet.eq.1.d0) then
         if (m.eq.k) djmk = 1.d0
         return
      endif

      if (cosbet.eq.-1.d0) then
         if (m.eq.is-k) djmk = dfloat(iv(j-m))
         return
      endif

      JMM = J-M
      JMK = J-K
      JPM = J+M
      JPK = J+K
      MPK = M+K

      IF (IS.EQ.1) THEN

         JPM = JPM-1
         JPK = JPK-1
         MPK = MPK-1

      ENDIF

      C2 = (1.d0+COSBET)/2
      S2 = (1.d0-COSBET)/2
      C  = SQRT(C2)
      S  = SQRT(S2)
      CS = C2/S2
      IA = MAX(0,-MPK)
      IE = MIN(JMM,JMK)
      X  = IV(JMM+IA) * C**(2*IA+MPK) * S**(JMM+JMK-2*IA) &
     &     *FI(JMM-IA)*FI(JMK-IA)*FI(MPK+IA)*FI(IA)
      do i = ia,ie
         DJMK = DJMK + X
         X = -X*CS*(JMM-I)*(JMK-I)/((I+1)*(MPK+I+1))
      enddo
      DJMK = DJMK*WF(JPM)*WF(JMM)*WF(JPK)*WF(JMK)
      return
!-end-djmk
      end

      double precision function wignei(j1,j2,j3,m1,m2,m3)

!=======================================================================
!
!     Calculates Wigner-coefficients for integer j- and m-values
!
!-----------------------------------------------------------------------
      USE VAPHFB_Par 
      implicit double precision (a-h,o-z)
!
      wignei=0.d0
      if (m1+m2+m3.ne.0) return
      i0 = j1+j2+j3+1
      if (igfv.lt.i0) stop 'in wignei: igfv too small'
      i1 = j1+j2-j3
      i2 = j1-m1
      i3 = j2+m2
      i4 = j3-j2+m1
      i5 = j3-j1-m2
      n2 = min0(i1,i2,i3)
      n1 = max0(0,-i4,-i5)
      if (n1.gt.n2) return
      do 11 n = n1,n2
   11 wignei = wignei + iv(n)*  &
     &         fi(n)*fi(i1-n)*fi(i2-n)*fi(i3-n)*fi(i4+n)*fi(i5+n)
      wignei = iv(i2+i3)*wignei*wfi(i0)*wf(i1)*wf(i2+i4)*wf(i3+i5)* &
     &         wf(i2)*wf(j1+m1)*wf(j2-m2)*wf(i3)*wf(j3-m3)*wf(j3+m3)
!
      return
!-end-WIGNEI
      end


        subroutine Num1_InfoMessage_LE20(value,Mess)
        implicit none
        character(len=20) :: Mess
        real*8 value
1       format(a20,'=',f7.3)
        write(*,*) '.....................'
        write(*,1)  Mess,value
        write(*,*) '.....................'

        end



        subroutine InfoMessage_LE20(Mess)
        implicit none
        character(len=20) :: Mess
1       format('..............:',a20)
        write(*,1) Mess

        end




      Double precision function dwignerI_gen(ttt,I,M,MP,beta)
!+---------------------------------------------------------------------+
!    This definition is different from the djmk
!    of Edmonds by a factor of (-1)^(M-MP)
!+---------------------------------------------------------------------+
!|   Wigner rotation matrix    (integer or half-integer                |
!|                                 angular momentum)                   |
!|                                                                     |
!|    I                 J-M'                            1/2            |
!|   d   , (beta) = (-1)    [(J+M)!(J-M)!(J+M')!(J-M')!]     ·         |
!|    M M                                                              |
!|          ---                  M+M'+2k              2J-M-M'-2k       |
!|          \    k  (cos(beta/2))        (sin(beta/2))                 |
!|         ·|(-1)  ------------------------------------------------    |
!|          /              k!(J-M-k)!(J-M'-k)!(M+M'+k)!                |
!|          ---                                                        |
!|           k                                                         |
!|                                                                     |
!|  k runs over all integer values for which the factorial arguments   |
!|  are non-negative                                                   |
!|                                                                     |
!|  IF 't'='i' ----> Integer angular momentum; I,M,MP are the actual   |
!|                                 values.                             |
!|                                                                     |
!|  IF 't'='h' ----> half-integer angular momentum; The actual values  |
!|                   for I,M,MP are I/2,M/2,MP/2                       |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|  Racah formula from Varshalovich                                    |
!|  (Quantum theory of angular momentum)                               |
!|                                                                     |
!+---------------------------------------------------------------------+

      Implicit real*8(A-H,O-Z),logical (L)
      Parameter (NFACM=150)
      character*1 ttt
      Dimension fact(0:NFACM)

      Save notcall,fact

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Factorials
      if(notcall.ne.61060) then
         fact(0) = 0.0d+00
         fact(1) = 0.0d+00
         do ifac =2,nfacm
            fact (ifac) = dlog(dfloat(ifac))  + fact (ifac-1)
         end do
         notcall = 61060
      end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      l1 = I.ge.IABS(M)
      l2 = I.ge.IABS(MP)

      lsr = l1.and.l2

      IF (ttt.eq.'i') THEN
       IplusM=I+M
       IminusM=I-M
       IplusMP=I+MP
       IminusMP=I-MP
       MplusMP=M+MP
       I2=2*I
      ELSE IF (ttt.eq.'h') THEN
       IplusM=(I+M)/2
       IminusM=(I-M)/2
       IplusMP=(I+MP)/2
       IminusMP=(I-MP)/2
       MplusMP=(M+MP)/2
       I2=I
      END IF

      if (lsr) then

      pre = 0.5d+00*(fact(IplusM)+fact(IminusM)+fact(IplusMP)+fact(IminusMP))

      pre = dexp (pre)

      phas = dfloat(1-2*mod(Iabs(IminusMP),2))

      dc   = dcos(0.5d+00*beta)
      ds   = dsin(0.5d+00*beta)

      fc   = dsign(1.0d+00,dc)
      fs   = dsign(1.0d+00,ds)

      if (fc*fs.lt.0.0d+00) then
         phas = phas*dfloat(1-2*mod(Iabs(MplusMP),2))
      end if

      dc   = max(dabs(dcos(0.5d+00*beta)),1.0d-20)  ! To prevent 0**0
      ds   = max(dabs(dsin(0.5d+00*beta)),1.0d-20)  !

      dlc  = dlog(dc)
      dls  = dlog(ds)

      imax = min(IminusM,IminusMP)
      imin = max(0,-(MplusMP))

      sum = 0.0d+00

      do k=imin,imax

         ss = fact(k)+fact(IminusM-k)+fact(IminusMP-k)+fact(MplusMP+k)
         ss = dfloat(MplusMP+2*k)*dlc+dfloat(I2-(MplusMP)-2*k)*dls - SS
         ss = dexp (ss)

         if(mod(k,2).eq.0) then
           sum = sum + ss
         else
           sum = sum - ss
         end if
        
      end do

      dwignerI_gen = phas*pre*sum
        
      else
      dwignerI_gen = 0.0d+00
      end if

      return
      end

      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

      SUBROUTINE ZCONST
!     ..........................................
      USE VAPHFB_PAR
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      COMMON IOUT,IREAD,IPRI
      COMMON /CHAT/ ZHAT(0:200)
      COMMON /CLOG/ ZLOG(100),ZGAMM2(100)
      COMMON /CONST/ ZAC,ZACM,ZPI,ZEPS
      ZAC=SQRT(2.0D0)
      ZACM=1.0D0/ZAC
      ZPI=4.0D0*DATAN(1.0D0)
      ZEPS=1.0D-12
      ZLOG(1)=0.0D0
      ZLOG(2)=0.0D0
      WN=1.0D0
      DO 1 I=3,100
      WN=WN+1.0D0
1     ZLOG(I)=ZLOG(I-1)+DLOG(WN)
      zgamm2(1)=dsqrt(zpi)
      zgamm2(2)=1.0d0
      do jj=3,100
         zgamm2(jj)=(jj-2)*zgamm2(jj-2)/2.0d0
      end do
      DO 3 S=0,200
3     ZHAT(S)=DSQRT(DBLE(S+1))
      RETURN
      END

!=======================================================================
      subroutine gfv 
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
      do m=1,igfvbc
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


       subroutine Real2Complex(akin,zakin,NLEV) 
       USE VAPHFB_PAR
       implicit none
       integer ii,jj,NLEV
       real*8     akin(NLEV,NLEV) 
       complex*16 zakin(NLEV,NLEV)

!        zone = (one,zero)
        do ii=1,NLEV
        do jj=1,NLEV
           zakin(jj,ii) = akin(jj,ii)*zone
        enddo
        enddo

        end

        integer function itconv(tz)
        real*8 tz
        itconv = 0
        if(abs(tz-1.d0) .lt. 1.d-3) itconv = 1    ! n
        if(abs(tz+1.d0) .lt. 1.d-3) itconv = -1   ! p
        if(itconv.eq.0) stop 'error in function itconv !'
        return
        end

        integer function mit(mt)
        integer mt

        if(mt .eq. -1) mit = 1   ! p
        if(mt .eq.  1) mit = 0   ! n
        if(abs(mt).ne.1) stop 'mt should be +1 or -1'

        return 
        end

        subroutine print_zmarray(LA1,LA2,LB1,LB2,LC1,LC2,LD1,LD2,ZE,ZA) 
        USE VAPHFB_PAR
        implicit none
        integer LA1,LB1,LC1,LD1
        integer LA2,LB2,LC2,LD2
        integer LA,LB,LC,LD
        complex*16 ZE,ZA
        real*8     AE 
        DIMENSION ZA(LA1:LA2,LB1:LB2,LC1:LC2,LD1:LD2)  ! it, lj, n1,n2
!        real*8, parameter :: CHOP=1.d-20
!        real*8  eps
!        eps=3.d-14
        if (abs(ZE) .le. CHOP) stop 'ZE is too small'
        do LA=LA1,LA2 
        do LB=LB1,LB2 
        do LC=LC1,LC2 
        do LD=LD1,LD2 
           AE= dreal(ZA(LA,LB,LC,LD)/ZE)
           if(abs(AE) .gt. CHOP) write(120,'(4i5,f12.8)') LA,LB,LC,LD,AE

        enddo
        enddo
        enddo
        enddo
        end
        
        integer function  SPB_twoj_lj(lj)
        implicit none
        integer lj
        SPB_twoj_lj = 1+ 2*int(lj/2) 
        end

        integer function SPB_l_lj(lj)
        implicit none
        integer lj
        SPB_l_lj = (lj+1)/2    
        end



      subroutine progress_bar( step  )
      implicit none
      integer :: step

      if ( step .ne. 0) then
     ! hold the backspace key down
       write(6,'(A)',advance='no') char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)&
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)
      flush 6
      end if

      write(6,'((A19),(I5))',advance='no') '    Completed(%):',min(100,step)
      flush 6

      end subroutine


        subroutine progress(j)
        implicit none
        integer(kind=4)::j,k
        character(len=17):: bar="???% |          |"
        write(unit=bar(1:3),fmt="(i3)") 10*j
        do k=1, j
           bar(6+k:6+k)="*"
        enddo
! print the progress bar.
        write(unit=6,fmt="(a1,a1,a17)") '+',char(13), bar
        return
        end subroutine progress


      subroutine zordls(n,ze)
!=======================================================================
!     orders a set of complex numbers according to their absolute value
!     from large to small values
!-----------------------------------------------------------------------
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      dimension ze(n)
!
      do 10 i = 1,n
         k  = i
         zp  = ze(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (abs(ze(j)).gt.abs(zp)) then
                  k = j
                  zp = ze(j)
               endif
   20       continue
            if (k.ne.i) then
               ze(k)  = ze(i)
               ze(i)  = zp
            endif
         endif
   10 continue
!
      return
!-end-ZORDLS
       end

!    ......................
      subroutine ordls(n,e)
!=======================================================================
!     orders a set of numbers according to their size
!     from large to small values
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension e(n)
!
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).gt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
            endif
         endif
   10 continue
!
      return
!-end-ORDLS
      end
