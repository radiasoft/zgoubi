      implicit double precision (a-h,o-z)

      character drctry(200)*132,zgDatFile(200)*132, txt132*132
      character cmmnd*300, fname*50
      integer debstr, finstr
      parameter(lunR=7,lunW=8)

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'NOW RUNNING PGM GENETEXLOG... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

C Create ./Log/log.tex and write head section
      cmmnd = 
     >'~/zgoubi/struct/tools/spin/xing_geneTexLog/geneTexIni'
      write(*,*) ' '
      write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
      call system(cmmnd)

C Get directory names
      open(unit=lunR,file='scanSpinResonances.Out2',err=99)
      read(lunR,*) txt132
      write(*,*) txt132
      j = 1
 1    continue
        read(lunR,*,end=10,err=10) i,drctry(i),zgDatFile(i)
        fname = zgDatFile(i)(debstr(zgDatFile(i)):finstr(zgDatFile(i)))
        write(*,*) i,drctry(i)(debstr(drctry(i)):finstr(drctry(i))),
     >         '   ',fname
        j = j + 1
        goto 1
 10   continue
      ifile = j-1

C Create a page for each xing
      do i = 1, ifile
        fname = './'//drctry(i)(debstr(drctry(i)):finstr(drctry(i)))//
     >    '/geneTexLog.out'
        open(unit=lunW,file=fname(debstr(fname):finstr(fname)),err=99)
        write(lunW,*) drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
        close(lunW)
        fname = drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
        cmmnd = 'cd '//fname(debstr(fname):finstr(fname))//' ; '//
     >  '~/zgoubi/struct/tools/spin/xing_geneTexLog/geneTexPage'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = '(cd '//fname(debstr(fname):finstr(fname))//' ; '//
     >  'cat logPage.tex >> ../Log/log.tex)'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
      enddo

C Rename all dirName/gnuplot_XXX.eps files to dirName/dirName_XXX.eps
C and cp to ./Log
      do i = 1, ifile
        fname = drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xingFull.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xingFull.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xingInit.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xingInit.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xingFinal.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xingFinal.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xing-xxp.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xing-xxp.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xing-yyp.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xing-yyp.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'mv '//fname(debstr(fname):finstr(fname))//
     >  '/gnuplot_xing-ldp.eps ' 
     >  //fname(debstr(fname):finstr(fname))//'/'
     >  //fname(debstr(fname):finstr(fname))//'_xing-ldp.eps'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
        cmmnd = 'cp '//fname(debstr(fname):finstr(fname))
     >  //'/*.eps ' 
     >  //'  ./Log'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)
      enddo

C End section of geneTexlog.tex
      cmmnd = 
     >'~/zgoubi/struct/tools/spin/xing_geneTexLog/geneTexEnd'
      write(*,*) ' '
      write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
      call system(cmmnd)

C 
        fname = drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
        cmmnd = '(cd ./Log ;  latex log ; xdvi log &)'
        write(*,*) ' '
        write(*,*) 'Pgm geneTexLog, execute : ',cmmnd
        call system(cmmnd)

      stop

 99   stop 'Pgm geneTexLog. Error open lunR'
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


      
      
