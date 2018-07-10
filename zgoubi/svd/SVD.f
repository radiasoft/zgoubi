      SUBROUTINE SVD(A,U,S,V,M,N)
      DOUBLE PRECISION A(M,N),U(M,M),VT(N,N),S(N),V(N,N)
!
! Program computes the matrix singular value decomposition.
! Using Lapack library.
!
! Programmed by sukhbinder Singh
! 14th January 2011
!


      DOUBLE PRECISION,ALLOCATABLE :: WORK(:)
      INTEGER LDA,M,N,LWORK,LDVT,INFO
      CHARACTER  JOBU, JOBVT

      JOBU='A'
      JOBVT='A'
      LDA=M
      LDU=M
      LDVT=N

      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))

      ALLOCATE(work(lwork))

      CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )

       DO I=1,2
         DO J=1,2
           V(J,I)=VT(I,J)
         END DO
       END DO
      END
