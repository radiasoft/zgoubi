! 1/ Create an ensemble of zgoubi_xxx.dat files for Froissard-Stora resonance xing, saved in respective xxx directories. 
! Starting data is calcStrength.out (contains list of resonances and their strengths),
! as obtained from prior executing /mad/tools/spin/calcStrength/calcStrength. 
! A zgoubi_xxx.dat is then built executing ~/zgoubi/toolbox/spin/xing_geneZgDat/geneZGDat4Xing_fromCalcStrength
! 2/ Run each zgoubi_xxx.dat in its own directory xxx. 
! 3/ executes various programs that compute resonance strengths from Sz vs. turn# etc.
! 4/ Execute various gnuplot.cmd types of files that produce graphs for monitoring, 
! 5/ Gather results from 3/ and 4/ into log.tex file

      implicit double precision (a-h,o-z)
      
      character radic*36, cmmnd*200, txtM*4, txtQz*15
      character txt300*300
      character zgDatFile(200)*50, drctry(200)*15
      integer debstr, finstr
      parameter(lunR=7,lunW=8,lunO=10,luntmp=14)

      character(100) dirTools
      parameter (dirTools=
     >   '~/zgoubi/toolbox/spin/resonanceXing/')

C scanSpinResonances.In is a copy of calcStrength.out 
C      call system('cp -f calcStrength.out scanSpinResonances.In')

      open(unit=lunR,file='scanSpinResonances.In',err=91)
      open(unit=lunO,file='scanSpinResonances.Out')

      write(lunO,fmt='(a)') 
     >   ' Data files created from calcStrength.out : '

C Generate zgoubi_xxx.dat file series
      ifile = 1
 1    continue
        read(lunR,fmt='(a)',end=18,err=99) txt300
C        close(lunR)
        open(unit=lunW,file='geneZGDat4Xing.In')
        write(lunW,*) txt300
        close(lunW)     

        read(txt300,*,end=18,err=99) 
     >    txtM, txtQz, aJn, aNn, znu, zmax, vkick
        ltxtM = finstr(txtM) - debstr(txtM) +1
        ltxtQz = finstr(txtQz) - debstr(txtQz) +1

        if(txtQz(debstr(txtQz):debstr(txtQz)+2) .eq. '|na') txtQz='000'

        drctry(ifile) = 'M'//txtM(debstr(txtM):debstr(txtM)+ltxtM-1)//
     >         'Qz'//txtQz(debstr(txtQz):debstr(txtQz)+ltxtQz-1)
        zgDatFile(ifile) = 'zgoubi_'
     >    //drctry(ifile)(debstr(drctry(ifile)):finstr(drctry(ifile)))
     >                           //'.dat'

C Temporary storage, for use by geneZGDat4Xing_fromCalcStrength   
        open(unit=luntmp,file='scanSpinResonances.tmp')
        write(luntmp,fmt='(i4,2(1X,a))') 
     >  ifile,drctry(ifile)(debstr(drctry(ifile)):finstr(drctry(ifile)))
     >  ,zgDatFile(ifile)(
     >     debstr(zgDatFile(ifile)):finstr(zgDatFile(ifile)))
        write(luntmp,fmt='(a4,1x,a15,1p,5(1x,e16.8),2a)') 
     >  txtM, txtQz, aJn, aNn, znu, zmax, vkick, 
     >  '  txtM, txtQz, aJn, aNn, znu, zmax, vkick ',
     >  ' cc''ed from geneZGDat4Xing.In'
        close(luntmp)

C Generate zgoubi_geneZGDat4Xing-Out.dat 
      cmmnd = dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_geneZgDat/geneZGDat4Xing_fromCalcStrength'
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)
         
C Create dedicated directory and 
C   moves zgoubi_geneZGDat4Xing-Out.dat there under dedicated name
        cmmnd = 'mkdir -p '//drctry(ifile)
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)
        cmmnd = 'cp -f zgoubi_geneZGDat4Xing-Out.dat '//
     >   zgDatFile(ifile)(
     >     debstr(zgDatFile(ifile)):finstr(zgDatFile(ifile)))
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        cmmnd = 'mv '//zgDatFile(ifile)//' '//drctry(ifile)
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)
        cmmnd = 'cp -f  scanSpinResonances.tmp  '//
     >  drctry(ifile)(
     >     debstr(drctry(ifile)):finstr(drctry(ifile)))
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)
        ifile = ifile+1        

        goto 1

 18   continue
      close(lunR)
      ifile = ifile -1
      write(*,*) ' ----------------------------------------------------'
      write(*,*) ' Reached the end of scanSpinResonances.In, ',
     > 'number of zgoubi.dat files created : ',ifile
      write(*,*) ' Names have been stored in scanSpinResonances.Out'
      write(*,*) 

C Save dirctory and file names for further use of xing_dataTreatment etc.
      do i = 1, ifile
        write(lunO,*) i,' ',drctry(i),' ',zgDatFile(i)
      enddo
      write(*,*) ' Now running this series of zgoubi.dat files... '
      write(*,*) ' ----------------------------------------------------'

C Run first created zgoubi_xxx.dat file. Keyword 'SYSTEM' at end of zgoubi_xxx.dat file will 
C launch run of the second one, and so forth. 
      open(unit=16,file='scanSpinResonances.count')
      read(16,*) i
      close(16) 
       cmmnd = 'cd '//
     >drctry(i)(
     >   debstr(drctry(i)):finstr(drctry(i)))
     >//' ; '
     >//'cp -f '//
     >zgDatFile(i)(
     >   debstr(zgDatFile(i)):finstr(zgDatFile(i)))
     >//'  zgoubi.dat'
C        write(*,*) ' Pgm scanSpinResonances, execute : ',cmmnd
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
      call system(cmmnd)
      cmmnd = 'cd '//
     >drctry(i)(
     >   debstr(drctry(i)):finstr(drctry(i)))
     >//' ; sed -i ''s/Data generated by searchCO/'
     >//drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
     >//'/g''   zgoubi.dat'
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
      call system(cmmnd)
      cmmnd = 'cd '//
     >drctry(i)(
     >   debstr(drctry(i)):finstr(drctry(i)))
     >//' ; ~/zgoubi/current/zgoubi/zgoubi < zgoubi.dat > echo'
      write(*,*) ' Pgm scanSpinResonances, now doing : ',
     >cmmnd(debstr(cmmnd):finstr(cmmnd))
        call system(cmmnd)

      stop

 99   continue
      ierr = 1
      write(*,*) ' '
      write(*,*) ' *** error during read in data file'
      write(*,*) ' '
      stop

 91   continue
      write(*,*) ' '
      write(*,*) 'PGM scanSpinResonances.f. '
     >  ,'Can''t open  scanSpinResonances.In'
      write(*,*) ' '
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

