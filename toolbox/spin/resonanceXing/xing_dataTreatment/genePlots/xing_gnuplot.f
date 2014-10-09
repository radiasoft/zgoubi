      character cmmnd*200, radic*90
      integer debstr, finstr
      character(100) dirTools
      parameter (dirTools=
     >   '/home/fmeot/zgoubi/toolbox/spin/resonanceXing/')

      cmmnd = '~/zgoubi/toolbox/fromBFai2Fai/fromBFai2Fai'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)

      radic = 'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_gnuplots'
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_xingFull.cmd'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_xingInit.cmd'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_xingFinal.cmd'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_xxp.cmd'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_yyp.cmd'
      write(*,*) ' command : ',cmmnd
      call system(cmmnd)
      cmmnd = radic(debstr(radic):finstr(radic))//
     >  '/gnuplot_ldp.cmd'
      write(*,*) ' command : ',cmmnd
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











