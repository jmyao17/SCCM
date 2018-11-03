!_______________________________________________________________________
      SUBROUTINE HNDIAG(NMAX,N,M,HH,NN,E,EN,FF,DD,GG,RR,WW,Z,NOS,IS,IFL)
!
!......................................................................C 
!     diagonalizes    HH*FF = E * NN*FF                                C
!     diagonalizes    NN*DD = EN * DD                                  C
!     eigenvalus of the norm EN<EPS are neglected                      C
!     IS = 0   only eigenvalues E,EN  and eigenvectors GG,DD           C
!     IS = 1   also eigenvectors  FF (possibly redundant)              C
!     IS = 2   also covariant components  WW = NN*FF                   C
!     IS = 3   also orthogonal wave functions RR = DD * GG             C
!              NN**(-1/2) * HH * NN**(-1/2) * RR  =  E  RR             C
!                                                                      C
!     HH(N,N),GG(N,N),NN(N,N),DD(N,N),FF(N,N),WW(N,N),RR(N,N)          C      
!     E(N)  eigenvalues of the generalized eigenvalue problem          C      
!     EN(N) eigenvalues of the norm matrix NN                          C
!     Z(N)  auxiliary field                                            C       
!                                                                      C
!......................................................................C
      USE VAPHFB_Par
      implicit real*8 (a-h,o-z)
      real*8 NN
!                                                                             
      dimension HH(NMAX,N),NN(NMAX,N),GG(NMAX,NMAX),DD(NMAX,N),FF(NMAX,N)
      DIMENSION WW(NMAX,N),RR(NMAX,N)
      dimension E(N),E1(N),EN(N),Z(N)
!------------------------------------                                                      
      EPS = GCM%eps !1.d-6 
      IFL = 0  
      do K=1,N
         do I=1,N
            DD(I,K)=-NN(I,K)
        enddo
      enddo
!-----diagonalization of NN
      CALL SDIAG(NMAX,N,DD,EN,DD,Z,1)
!      
      do I=1,N
         EN(I)=-EN(I)
      enddo
!     ............... number of natural state is defined by the cutoff in the eigen values of norm 
!--- find out all the eigenvalues which are large enough compared with the largest one,
!    and drop the small ones
      do 1 I=1,N
         IF (EN(I)/EN(1).LT.EPS) GOTO 2
    1    CONTINUE
    2 M=I-1
!     ............... number of natural state is given by NOS
      M=min(M,NOS)
      IF (M.LE.0) THEN
         WRITE(6,*) ' SUBR. HNDIAG: NO EIGENVALUE OF N LARGER THAN EPS'
         IFL=1
         RETURN
      ENDIF
!----- up here, all left eigenvalues are large enough, the number is M
!      parepare the new Hamiltonian in the collective subspace
!-----
!----- DD = u
      do K=1,M
         S=1./SQRT(EN(K))
         do I=1,N
           RR(I,K)=DD(I,K)*S
         enddo !I
      enddo !K
!-------------
      do K=1,M
         do I=1,N
            S=0.
            do L=1,N
               S=S+HH(I,L)*RR(L,K)
            enddo !L
            E1(I)=S
         enddo !I
         do J=K,M
            S=0.
            do L=1,N
               S=S+RR(L,J)*E1(L)
            enddo ! L
               GG(J,K)=S
         enddo !J 
       enddo !K
  100 format(5x,12f10.5)

!---- diagonalization of the new Hamiltonian
! 
      CALL sdiag(NMAX,M,GG,E,GG,Z,1) ! dim(E) = M 
!
!     calculation of the contra-variant components FF
      IF (IS.EQ.0) RETURN
      CALL MAB(NMAX,NMAX,NMAX,N,M,M,RR,GG,FF,1,0)
!
!     calculation of the co-variant components WW (overlaps)
      IF (IS.LT.2) RETURN
      CALL MAB(NMAX,NMAX,NMAX,N,N,M,NN,FF,WW,1,0)
!
!     BESTIMMUNG DER ORTHOGONALEN WELLENFUNKTIONEN RR
      IF (IS.LT.3) RETURN
!     calculation of the orthogonal wave functions RR
      CALL MAB(NMAX,NMAX,NMAX,N,M,M,DD,GG,RR,1,0)
!-----
      RETURN
!-end-HNDIAG
      END
!_______________________________________________________________________
      subroutine mab(ma,mb,mc,n1,n2,n3,a,b,c,iph,is)
!
!......................................................................
!
      implicit real*8 (a-h,o-z)
!
      dimension a(ma,ma),b(mb,mb),c(mc,mc)
!      dimension a(n1,n2),b(n2,n3),c(n1,n3)
!
      if (is.eq.0) then
         if (iph.lt.0) then
            do i = 1,n1
            do k = 1,n3
               s = 0.0d0
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k) = - s
            enddo
            enddo
         else
            do i = 1,n1
            do k = 1,n3
               s = 0.0d0
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k) = s
            enddo
            enddo
         endif
      else
         if (iph.lt.0) then
            do i = 1,n1
            do k = 1,n3
               s = c(i,k)
               do l = 1,n2
                  s = s-a(i,l)*b(l,k)
               enddo
               c(i,k) = s
            enddo
            enddo
        else
            do i = 1,n1
            do k = 1,n3
               s = c(i,k)
               do l = 1,n2
                  s = s+a(i,l)*b(l,k)
               enddo
               c(i,k)=s
            enddo
            enddo
         endif
      endif
      return
!-end-MAB
      end

