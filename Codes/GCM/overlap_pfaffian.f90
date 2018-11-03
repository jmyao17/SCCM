   SUBROUTINE Overlap_Pfaffian(dim, U0mat, V0mat, U1mat, V1mat, Rotmat, overlap, STAT)
   USE vaphfb_par 
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: dim
   COMPLEX(DP), DIMENSION(1:dim,1:dim), INTENT(IN) :: U0mat, V0mat, U1mat, V1mat, Rotmat
   COMPLEX(DP), INTENT(OUT) :: overlap
   INTEGER, INTENT(OUT) :: STAT

   REAL(DP), DIMENSION(1:dim/2) :: u0_can, u1_can, v0_can, v1_can
   COMPLEX(DP), DIMENSION(1:dim,1:dim) :: D0mat, C0mat, D1mat, C1mat   
   COMPLEX(DP), DIMENSION(1:2*dim,1:2*dim) :: Wmat, WORK, Mmat

   REAL(DP) :: Norm0, Norm1
   COMPLEX(DP) :: Delta, detD0, detC0, detD1, detC1, Constant_W, Constant_M, overlap_W, overlap_M
   
   COMPLEX(DP), DIMENSION(1:dim,1:dim) :: V0matinv, V1matinv, UV0mat, UV1mat, UV0bar, UV1bar, Vnew
   INTEGER :: i, IH, STAT0, STAT1, j
   COMPLEX(DP) :: Pfaffian_W, Pfaffian_M
   ! compute the overlap  < U0V0 | U1V1 >
   
   ! dimension check
   IF( MOD(dim, 2) .NE. 0) THEN
      STAT = -1
      STOP "overlap_pfaffian: dim should be an even number"
   END IF
        
   ! bloch-messiah
   PRINT *, "Bloch Messiah for ini"
   CALL BlochMessiah(dim,U0mat,V0mat,u0_can,v0_can,D0mat,C0mat,Vnew,STAT0)
   PRINT *, "Bloch Messiah for fin"   
   CALL BlochMessiah(dim,U1mat,V1mat,u1_can,v1_can,D1mat,C1mat,Vnew,STAT1)
   IF(STAT0 .NE. 0 .OR. STAT1 .NE. 0) THEN
      WRITE(0,*) "overlap_pfaffian: Bloch Messiah failed"
      STAT = -1
      RETURN
   END IF

   PRINT *, "Determinant"
   CALL COMPLEX_DetofMat(dim, D0mat, detD0)
   CALL COMPLEX_DetofMat(dim, C0mat, detC0)
   CALL COMPLEX_DetofMat(dim, D1mat, detD1)
   CALL COMPLEX_DetofMat(dim, C1mat, detC1)
   PRINT *, "determinant done"
   PRINT *, "detD0 = ", detD0
   PRINT *, "detC0 = ", detC0
   PRINT *, "detD1 = ", detD1
   PRINT *, "detC1 = ", detC1
   
   Norm0 = 0.0D0; Norm1 = 0.0D0
   DO i = 1, dim/2
      Norm0 = Norm0 + LOG(v0_can(i))
      Norm1 = Norm1 + LOG(v1_can(i))
      PRINT *, "i  = ", i, v0_can(i) , v1_can(i)      
   END DO
   PRINT *, "Norm0 = ", EXP(Norm0)
   PRINT *, "Norm1 = ", EXP(Norm1)
   
   Delta = CONJG(detD0) * detD1 * detC0 * CONJG(detC1)
   PRINT *, "Delta", Delta
   Constant_W = (-1)**(dim/4) * Delta * EXP(Norm0 + Norm1)
   PRINT *, "Constant", Constant_W
   Constant_M = (-1)**(dim/4) * EXP(-(Norm0 + Norm1))
   ! UV^{-1}
   PRINT *, "complex inverse"
   CALL COMPLEX_INVERSE(dim, V0mat, V0matinv, STAT)
   PRINT *, "STAT = ", STAT
   CALL COMPLEX_INVERSE(dim, V1mat, V1matinv, STAT)
   PRINT *, "STAT = ", STAT
!   UV0mat(1:dim,1:dim) = MATMUL(U0mat(1:dim,1:dim), V0matinv(1:dim,1:dim))
!   UV1mat(1:dim,1:dim) = MATMUL(U1mat(1:dim,1:dim), V1matinv(1:dim,1:dim))

   UV0bar(1:dim,1:dim) = 0.0D0
   UV1bar(1:dim,1:dim) = 0.0D0
   DO i = 1, dim/2
      UV0bar(2*i-1,2*i) = - u0_can(i) / v0_can(i) 
      UV0bar(2*i,2*i-1) =   u0_can(i) / v0_can(i)       
   END DO
   PRINT *, "Skewcheck UV0bar"
   CALL SkewCheck(dim, UV0bar, STAT)
   IF(STAT .NE. 0) WRITE(0,*) "UV0bar not skew"
      DO i = 1, dim/2
      UV1bar(2*i-1,2*i) = - u1_can(i) / v1_can(i) 
      UV1bar(2*i,2*i-1) =   u1_can(i) / v1_can(i)       
   END DO
   PRINT *, "Skewcheck UV1bar"
   CALL SkewCheck(dim, UV1bar, STAT)
   IF(STAT .NE. 0) WRITE(0,*) "UV1bar not skew"

   UV0mat(1:dim,1:dim) = MATMUL( D0mat(1:dim,1:dim), MATMUL( UV0bar(1:dim,1:dim), &
        & TRANSPOSE( D0mat(1:dim,1:dim))))
   UV1mat(1:dim,1:dim) = MATMUL( D1mat(1:dim,1:dim), MATMUL( UV1bar(1:dim,1:dim), &
        & TRANSPOSE( D1mat(1:dim,1:dim))))
   PRINT *, "Skewcheck UV0"
   CALL SkewCheck(dim, UV0mat, STAT)
   IF(STAT .NE. 0) WRITE(0,*) "UV0 not skew"
   PRINT *, "Skewcheck UV1"
   CALL SkewCheck(dim, UV1mat, STAT)
   IF(STAT .NE. 0) WRITE(0,*) "UV1 not skew"
   

   ! M matrix
   Mmat(1:2*dim,1:2*dim) = 0.0D0
   
   Mmat(    1:dim,      1:dim)   = MATMUL( &
        & TRANSPOSE(V0mat(1:dim,1:dim)), U0mat(1:dim,1:dim))
   Mmat(    1:dim,  1+dim:2*dim) = MATMUL ( &
        & TRANSPOSE(V0mat(1:dim,1:dim)), MATMUL(Rotmat(1:dim,1:dim), CONJG(V1mat(1:dim,1:dim))))
   Mmat(1+dim:2*dim,    1:dim) = - MATMUL ( &
        & CONJG(TRANSPOSE(V1mat(1:dim,1:dim))), MATMUL(TRANSPOSE(Rotmat(1:dim,1:dim)), V0mat(1:dim,1:dim)))
   Mmat(1+dim:2*dim,1+dim:2*dim) = MATMUL( &
        & CONJG(TRANSPOSE(U1mat(1:dim,1:dim))), CONJG(V1mat(1:dim,1:dim)))
   
   ! W matrix
   
   Wmat(1:2*dim,1:2*dim) = 0.0D0

   Wmat(    1:dim,      1:dim)   =  TRANSPOSE(CONJG(UV1mat(1:dim,1:dim)))
   Wmat(    1:dim,  1+dim:2*dim) = -TRANSPOSE(Rotmat(1:dim,1:dim))
   Wmat(1+dim:2*dim,    1:dim)   =  Rotmat(1:dim,1:dim)
   Wmat(1+dim:2*dim,1+dim:2*dim) =  UV0mat(1:dim,1:dim)
   PRINT *, "Skew check Wmat"
   CALL SkewCheck(2*dim, Wmat, STAT)
   IF(STAT .NE. 0) WRITE(0,*) "Wmat not skew"
   
   PRINT *, "dim = ", dim
   PRINT *, "U1 V1^{-1}"   
   DO i = 1, 10
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 1, 10)
   END DO
      PRINT *, ""
   DO i = 11, 20
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 11, 20)
   END DO
      PRINT *, ""
   DO i = 21, 30
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 21, 30)
   END DO
      PRINT *, ""
   DO i = 31, 40
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 31, 40)
   END DO
   PRINT *, "-Rotmat^{T}"
   DO i = 1, 10
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 1+dim, 10+dim)
   END DO
   PRINT *, "Rotmat"
   DO i = dim+1, 10+dim
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 1, 10)
   END DO
   PRINT *, "U0 V0^{-1}"
   DO i = 1+dim, 10+dim
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 1+dim, 10+dim)
   END DO
      PRINT *, ""
   DO i = 11+dim, 20+dim
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 11+dim, 20+dim)
   END DO
   PRINT *, ""
   DO i = 21+dim, 30+dim
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 21+dim, 30+dim)
   END DO
   PRINT *, ""
   DO i = 31+dim, 40+dim
      WRITE(*,'(10ES15.5)') (ABS(Wmat(i,j)), j = 31+dim, 40+dim)
   END DO
   ! check skew / construct the matrix from D and C of BM.
   Pfaffian_W = 0.0D0
   Pfaffian_M = 0.0D0
   IH = 2 ! or 2
   PRINT *, "ZPfaffianH" 
   CALL ZPfaffianH (Wmat,2*dim,2*dim,Pfaffian_W,WORK,IH)
   PRINT *, "PfaffianW = ", Pfaffian_W
   PRINT *, "ConstantW = ", Constant_W
   overlap_W = Constant_W * Pfaffian_W
   PRINT *, "overlap_W = ", overlap_W 
   CALL ZPfaffianH (Mmat,2*dim,2*dim,Pfaffian_M,WORK,IH)
   PRINT *, "PfaffianM = ", Pfaffian_M
   PRINT *, "ConstantM = ", Constant_M
   overlap_M = Constant_M * Pfaffian_M
   PRINT *, "overlap_M = ", overlap_M
   overlap = overlap_M
   
   STAT = 0
   
   RETURN
   
 END SUBROUTINE Overlap_Pfaffian
 
