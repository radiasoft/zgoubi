      FUNCTION ISNUM(S)
      IMPLICIT double precision (a-h,o-z)
      CHARACTER(*) S
      LOGICAL ISNUM
      INTEGER I, K, N
      integer debstr, finstr
      logical alrdy
      data alrdy / .false. /
      isnum = .false.
      N = LEN(s(debstr(s):finstr(s)))
      I = 1
      K = ICHAR(S(I:I))
      IF (K .lt. 48) THEN
        IF (.not.(K.eq.43 .or. K.eq.45)) THEN   ! + or -
          ISNUM = .FALSE.
          RETURN
        endif
      elseIF (K .gt. 57) THEN
        ISNUM = .FALSE.
        RETURN
      ENDIF
      DO I = 2, N
        K = ICHAR(S(I:I))
        IF (K .lt. 48) THEN
          IF (K.ne.46) THEN    !  a dot
            ISNUM = .FALSE.
            RETURN
          else
            if(alrdy) then
              ISNUM = .FALSE.
              RETURN
            else
              alrdy = .true.
            endif 
          endif
        elseIF (K .gt. 57) THEN
          ISNUM = .FALSE.
          RETURN
        ENDIF
      ENDDO
      ISNUM = .TRUE.
      RETURN
      END 
