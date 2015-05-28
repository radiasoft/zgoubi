      FUNCTION ISNUM(S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISNUM
      CHARACTER(*) S
      INTEGER I, K, N
      INTEGER DEBSTR, FINSTR
      LOGICAL ALRDY
      DATA ALRDY / .FALSE. /
      ISNUM = .FALSE.
      S = S(DEBSTR(S):FINSTR(S))
      N = LEN(S(DEBSTR(S):FINSTR(S)))
      I = 1
      K = ICHAR(S(I:I))
      IF (K .LT. 48) THEN
        IF (.NOT.(K.EQ.43 .OR. K.EQ.45)) THEN   ! + OR -
          ISNUM = .FALSE.
          GOTO 88
        ENDIF
      ELSEIF (K .GT. 57) THEN
        ISNUM = .FALSE.
        GOTO 88
      ENDIF
      DO I = 2, N
        K = ICHAR(S(I:I))
        IF (K .LT. 48) THEN
          IF (K.NE.46) THEN    !  A DOT
            ISNUM = .FALSE.
            GOTO 88
          ELSE
            IF(ALRDY) THEN
              ISNUM = .FALSE.
              GOTO 88
            ELSE
              ALRDY = .TRUE.
            ENDIF 
          ENDIF
        ELSEIF (K .GT. 57) THEN
          ISNUM = .FALSE.
          GOTO 88
        ENDIF
      ENDDO
      ISNUM = .TRUE.

 88   CONTINUE
      RETURN
      END 
