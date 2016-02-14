      character(200) txt, txt2
      dimension surf(200), pp0(200)
      logical strcon
      data y / 0. /

      open(unit=1,file='searchDA.out')
      open(unit=2,file='surface.out')
        
      do i = 1, 200
        surf(i) = 0.
      enddo
      smax = -1e10

 1    continue

        read(1,fmt='(a)',err=10,end=10) txt

        if(   (.not. strcon(txt,'#',is) )
     >  .and. (.not. strcon(txt,'dx-',is))  ) then

          read(1,fmt='(a)',err=10,end=10) txt2

          if(   (.not. strcon(txt2,'#',is) )
     >    .and. (.not. strcon(txt2,'dx-',is))  ) then

             read(txt,*,err=1,end=1) x1, y1, pp, xmom
             read(txt2,*,err=1,end=1) x2, y2, ppb, xmomb
             mom = nint(xmom)
             surf(mom) = surf(mom) + 0.5*(y1+y2)*(x2-x1)
             lastm = mom
             pp0(mom) = pp
             if (surf(mom) .gt. smax) smax = surf(mom)
             lmom = mom
             write(*,*) x1,x2,0.5*(y1+y2), pp, mom,surf(mom),smax,lastm
          endif
        endif

      goto 1

 10   continue

      do mom =1, lmom
         write(2,*) pp0(mom), surf(mom), surf(mom)/smax*pp0(mom),mom
      enddo        
      close(2)
      stop
      end

      FUNCTION STRCON(STR,STR2,
     >                         IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRCON
      CHARACTER(*) STR, STR2
C     ---------------------------------------------------------------
C     .TRUE. if the string STR contains the string STR2 at least once
C     IS = position of first occurence of STR2 in STR 
C     (i.e.,STR(IS:IS+LEN(STR2)-1)=STR2)
C     ---------------------------------------------------------------
      INTEGER DEBSTR,FINSTR
      LNG2 = LEN(STR2(DEBSTR(STR2):FINSTR(STR2)))
      IF(LEN(STR).LT.LNG2 .OR.
     >   (DEBSTR(STR).EQ.0 .AND. FINSTR(STR).EQ.0)
     >     ) GOTO 1
      DO I = DEBSTR(STR), FINSTR(STR)-LNG2+1
        IF( STR(I:I+LNG2-1) .EQ. STR2 ) THEN
          IS = I 
          STRCON = .TRUE.
          RETURN
        ENDIF
      ENDDO
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
      FUNCTION DEBSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER(*) STRING
C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------
      II=0
      LENGTH=LEN(STRING)
      IF(LENGTH.LE.0) GOTO 99
1     CONTINUE
        II=II+1
        IF (STRING(II:II) .EQ. ' ') THEN
          IF(II .GE. LENGTH) THEN
            II = 0
            GOTO 99
          ELSE
            GOTO 1
          ENDIF
        ENDIF
 99   CONTINUE
      IF(II .EQ. 0) THEN
        STRING = ' '
        II = 1
      ENDIF
      DEBSTR = II
      RETURN
      END
      FUNCTION FINSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FINSTR
      CHARACTER(*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) GOTO 99
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

 99   CONTINUE

      IF(FINSTR.EQ.0) THEN
        FINSTR = 1
        STRING = ' '
      ENDIF
      RETURN
      END
