      character cmmnd*132 
      character(100) dirTools
      parameter (dirTools=
     >   '/home/meot/zgoubi/struct/tools/spin/resonanceXing/')
      integer debstr, finstr

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'NOW RUNNING PGM GENEPLOTS... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

      cmmnd = 
     >'~/zgoubi/struct/tools/b_fai2Fai/fromBFai2Fai'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_xingFull.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_xingInit.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_xingFinal.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_xxp.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_yyp.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)
      
      cmmnd =  'gnuplot < '//
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_dataTreatment/genePlots/'//
     >'gnuplot_ldp.cmd'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
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

