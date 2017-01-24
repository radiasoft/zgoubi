      REAL*8 FUNCTION GETDET(A, N)
C A function written in FORTRAN77 to calculate determinant of a square matrix

C Passed parameters:
C A = the matrix
C N = dimension of the square matrix
       
C A modification of a code originally written by Ashwith J. Rego, available from http://www.dreamincode.net/code/snippet1273.htm
C Modified by Syeilendra Pramuditya, available from http://wp.me/p61TQ-zb
C Last modified on January 13, 2011 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ELEM(N,N),A(N,N)
      REAL*8 M, TEMP
      INTEGER I, J, K, L
      LOGICAL DETEXISTS

      DO I=1,N
      DO J=1,N
      ELEM(I,J)=A(I,J)
      END DO
      END DO
        
      DETEXISTS = .TRUE.
      L = 1
      !CONVERT TO UPPER TRIANGULAR FORM
      DO K = 1, N-1
	  IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
	  DETEXISTS = .FALSE.
	  DO I = K+1, N
	  IF (ELEM(I,K).NE.0.0) THEN
      
	  DO J = 1, N
	  TEMP = ELEM(I,J)
	  ELEM(I,J)= ELEM(K,J)
	  ELEM(K,J) = TEMP
	  END DO
      
	  DETEXISTS = .TRUE.
	  L=-L
	  EXIT
      
	  END IF
      
	  END DO
	  IF (DETEXISTS .EQV. .FALSE.) THEN
	  GETDET = 0
	  RETURN
	  END IF
      END IF
      DO J = K+1, N
	  M = ELEM(J,K)/ELEM(K,K)
	  DO I = K+1, N
	  ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
	  END DO
	  END DO
      END DO
	
      !CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
      GETDET = L
      DO I = 1, N
	  GETDET = GETDET * ELEM(I,I)
      END DO
	
      END        
