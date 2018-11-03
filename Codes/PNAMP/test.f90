         subroutine test()


        implicit none
        integer ii,jj
        integer INFO2,ndim
!        real*8 aa(3,3),ee(3),ez(3)
        real*8 aa(4,4),ee(4),ez(4),bb(4,4)

       aa = 0.d0

       ndim = 4
       aa(1,1) = 2
       aa(2,2) = 2
       aa(1,2) = 1
       aa(2,1) = 1

       aa(3,3) = 3
       aa(4,4) = 3
       aa(3,4) = 1
       aa(4,3) = 1

!  DSYEV: computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A        
!  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )            
!  ----------------------------
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!  A       (input/output) DOUBLE PRECISION array
!  W       (output) eigenvalues

       call dsyev('V','U',NDIM,aa,NDIM,ee,bb,3*NDIM-1,INFO2)
!       call sdiag(ndim,ndim,aa,ee,aa,ez,0)

       do ii=1,ndim
       do jj=1,ndim

          write(*,*) ii,jj,aa(jj,ii)

        enddo
          
        enddo
        do ii=1,ndim
          write(*,*) ii, ee(ii)
        enddo

        stop
       end

