      character(200) cmmnd
      character(100) dirTools
      parameter (dirTools=
     >   '/home/fmeot/zgoubi/toolbox/spin/resonanceXing/')
      integer debstr, finstr

C----------------------------
C Launch tracking procedure      
      open(unit=8,file='scanSpinResonances.count')
      i = 1
      write(8,fmt='(i4,a)') i,' = xing simulation number'
      close(8)
      cmmnd = dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_scanResonances/scanSpinResonances'
      write(*,fmt='(2a)') ' Pgm scanSpinResonances_launch, now doing: ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
      call system(cmmnd)

cC----------------------------
cC Creates log.tex 
c      cmmnd = '/home/fmeot/zgoubi/toolbox/spin/xing_geneTexLog/'
c     >  //'geneTexLog'
c      call system(cmmnd)      

C-----------------------------
C 
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' ////////////////////////////////////////////////////'
      write(*,*) ' PGM scanSpinResonances_launch. '
      write(*,*) ' xing_scanResonances/scanSpinResonances now launched,'
      write(*,*) ' "ps -u meot | grep zgoubi" allows checking that'
     > ,' the zgoubixxx.dat series is running... '
      write(*,*) ' ////////////////////////////////////////////////////'
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
