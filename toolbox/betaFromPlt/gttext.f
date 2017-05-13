      FUNCTION GTTEXT(NRES,LUNR,TXT,
     >                              TXTIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GTTEXT
      CHARACTER(*) TXT
      CHARACTER(*) TXTIN
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR

      READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN

      DOWHILE(.NOT.STRCON(TXTIN,TXT(DEBSTR(TXT):FINSTR(TXT)),
     >                                                       IS))
        READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN
      ENDDO

      GTTEXT = .TRUE.
      GOTO 10

 99   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - ERR upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - EOF upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 10   CONTINUE
      RETURN
      END
