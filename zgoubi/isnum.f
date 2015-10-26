      FUNCTION ISNUM(S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISNUM
      CHARACTER(*) S
      INTEGER I, K, N
      INTEGER DEBSTR, FINSTR
      LOGICAL ALRDY

      ALRDY = .FALSE. 

      ISNUM = .FALSE.
      S = S(DEBSTR(S):FINSTR(S))
      N = LEN(S(DEBSTR(S):FINSTR(S)))

      I = 1
      K = ICHAR(S(I:I))
      IF (K .LT. 48) THEN
        IF (.NOT.(K.EQ.43 .OR. K.EQ.45 .OR. K.EQ.46)) THEN   ! + OR - OR .
          ISNUM = .FALSE.
          GOTO 88
        ELSEIF(K.EQ.46) THEN
          ALRDY = .TRUE.
        ENDIF
      ELSEIF (K .GT. 57) THEN     ! 49-57 : 1-9
        ISNUM = .FALSE.
        GOTO 88
      ENDIF

      DO I = 2, N
        K = ICHAR(S(I:I))
        IF (K .LT. 48) THEN
          IF (K.NE.46) THEN    !  46 : DOT
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
          IF(  S(I:I).EQ.'e' .OR. S(I:I).EQ.'E' 
     >    .OR. S(I:I).EQ.'d' .OR. S(I:I).EQ.'D') THEN
            GOTO 77
          ELSE
            ISNUM = .FALSE.
            GOTO 88
          ENDIF
        ENDIF
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      ENDDO
      ISNUM = .TRUE.
      goto 88

 77   CONTINUE
      I = I+1
      K = ICHAR(S(I:I))
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      IF (K .LT. 48) THEN
        IF (.NOT.(K.EQ.43 .OR. K.EQ.45)) THEN   ! + OR - 
          ISNUM = .FALSE.
          GOTO 88
        ENDIF
      ELSEIF (K .GT. 57) THEN
        ISNUM = .FALSE.
        GOTO 88
      ENDIF
      J = I + 1
      DO I = J, N
        K = ICHAR(S(I:I))
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
        IF (K .LT. 48) THEN
            ISNUM = .FALSE.
            GOTO 88
        ELSEIF (K .GT. 57) THEN
            ISNUM = .FALSE.
            GOTO 88
        ENDIF
C      write(*,*) ' n, i, char, isnum ',n, i, ICHAR(S(I:I)),s(i:i)
      ENDDO
      ISNUM = .TRUE.


 88   CONTINUE

      RETURN
      END 
