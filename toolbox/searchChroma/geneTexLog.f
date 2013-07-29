      character*120 cmmnd

      character(30) tilde
      parameter (tilde = '/home/owl/fmeot')
      integer debstr, finstr

C Create ./Log if does not exist already
      cmmnd = 'mkdir -p Log_Chroma'
      call system(cmmnd)

      cmmnd = 'cp zgoubi_searchCO-*.dat Log_Chroma'
      call system(cmmnd)

      cmmnd = tilde(debstr(tilde):finstr(tilde))
     >  //'/zgoubi/struct/tools/searchChroma/geneTexIni'
      write(*,*) '------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)
      cmmnd = tilde(debstr(tilde):finstr(tilde))
     >  //'/zgoubi/struct/tools/searchChroma/geneTexBody'
      write(*,*) '-------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)
      cmmnd = tilde(debstr(tilde):finstr(tilde))
     >  //'/zgoubi/struct/tools/searchChroma/geneTexEnd'
      write(*,*) '------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      cmmnd = 
     >'cd Log_Chroma ; latex log_Chroma ; '
     >//'latex log_Chroma ; dvipdf log_Chroma'
      write(*,*) '-----------------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      stop
      end

      FUNCTION DEBSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)
C      LENGTH=LEN(STRING)+1
1     CONTINUE
        DEBSTR=DEBSTR+1
C        IF(DEBSTR .EQ. LENGTH) RETURN
C        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') GOTO 1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .EQ. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END
      FUNCTION FINSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER FINSTR
      CHARACTER * (*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) RETURN
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END

