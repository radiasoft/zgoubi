C dataTreatment pgm treats data in M*Qz* directories. 
C This assumes prior execution of scanSpinResonances_launch.f in M*Qz*s' parent directory.
C dataTreatment is to be launched from the M*Qz*s' parent directory, it will 
C find (and get !) M*Qz*s' names from scanSpinResonances.Out2 storage file. 

      implicit double precision (a-h,o-z)

      character txtM*4, txtQz*15, txt1 
      character drctry(200)*132,zgDatFile(200)*132
      character txt132*132,txt12*12,txt160*160
      parameter(lunR=7,lunW=8,lunMQ=9)

      character cmmnd*300, fname*50
      integer debstr, finstr
      logical strcon
      character typ2*12

      character(100) dirTools
      parameter (dirTools=
     >   '~/zgoubi/toolbox/spin/resonanceXing/')

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'NOW RUNNING PGM DATATREATMENT... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

C      open(unit=lunR,file='scanSpinResonances.Out2',err=99)
      open(unit=lunR,file='scanSpinResonances.Out',err=99)
      read(lunR,*) txt132
      write(*,*) txt132
      j = 1
 1    continue
        read(lunR,*,end=10,err=10) i,drctry(i),zgDatFile(i)
        write(*,*) i,drctry(i),zgDatFile(i)
        j = j + 1
        goto 1
 10   continue
      ifile = j-1

C Execute 'tunesFromFai' pgm :  computation of tunes and emittances 
      do i = 1, ifile
C Get number of turns, from zgoubi.res
        fname = 
     >  './'//drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
     >  //'/zgoubi.res'
        open(unit=34,file=fname,err=122)
 22     continue
          read(34,fmt='(a)',end=102) txt132
          txt132 = txt132(debstr(txt132):132)

          if(strcon(txt132,'''REBELOTE''',10,
     >                                     IS)) then
              read(34,*) npass
              goto 102
          endif
          goto 22
 102      continue
        close(34)
        fname = 
     >  './'//drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
     >  //'/tunesFromFai.In'
        open(unit=34,file=fname)
        write(34,*) int(npass/2*.8d0), int(npass/2*1.2d0)
     >                                      ,'  ! turn# range'
        write(34,*) '0.51, 1., 0.51, 1., 0., 1. '
     >                                  //'  ! tune boundaries'
        write(34,*) ' 2000, 2000, 2000, '
     >                                  //' !  spectrum sampling'
        write(34,*) ' n '//'  ! do not store spectra '
        write(34,*) ' y '//'  ! compute tunes '
        close(34)
        cmmnd = 'cd '//
     >  drctry(i)(debstr(drctry(i)):finstr(drctry(i)))//' ; '//
     >  ' ~/zgoubi/toolbox/tunesFromFai/tunesFromFai'
        write(*,*) ' '
        write(*,*) 'Pgm dataTreatment, execute : ',cmmnd
        call system(cmmnd)
      enddo

 122  continue

C Execute avrgSzFromFai pgm :  computation of average Sz over first and over last turns
      do i = 1, ifile
        cmmnd = 'cd '//
     >  drctry(i)(debstr(drctry(i)):finstr(drctry(i)))//' ; '//
     >  ' pwd ; '//
     >  dirTools(debstr(dirTools):finstr(dirTools))
     >  //'xing_dataTreatment/avrgSzFromFai/avrgSzFromFai'
        write(*,*) ' '
        write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
        call system(cmmnd)
      enddo

C Execute computeStrength pgm
      do i = 1, ifile
        cmmnd = 'cd '//
     >  drctry(i)(debstr(drctry(i)):finstr(drctry(i)))//' ; '//
     >  ' pwd ; '//
     >  dirTools(debstr(dirTools):finstr(dirTools))
     >  //'xing_dataTreatment/computeStrength/computeStrength'
        write(*,*) ' '
        write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
        call system(cmmnd)
      enddo

C Makes presen-table output
      open(unit=lunW,file='dataTreatment.out')
      open(unit=lunMQ,file='scanSpinResonances.In')

      do i = 1, ifile
        fname = 
     >  './'//drctry(i)(debstr(drctry(i)):finstr(drctry(i)))
     >  //'/computeStrengths.out'

        open(unit=34,file=fname,err=123)
        read(34,*) txt1, typ2
        read(34,fmt='(a)') txt160
        close(34)

        if    (typ2.eq.'intrinsic') then
          read(lunMQ,*) txtM, txtQz, dum, strMAD
          txt160 = txtM(debstr(txtM):finstr(txtM))//
     >    txtQz(debstr(txtQz):debstr(txtQz))//'Qz      '//
     >    txt160(debstr(txt160):finstr(txt160))
          write(txt12,fmt='(1p,g12.4)') strMAD
          txt160 = 
     >    txt160(debstr(txt160):finstr(txt160))
     >    //'  &  '//txt12//'  \\'
        elseif(typ2.eq.'imperfection') then
          read(lunMQ,*) txtM, txtQz, strMAD, dum, dum, zco
          txt160 = txtM(debstr(txtM):finstr(txtM))//
     >    '           '//
     >    txt160(debstr(txt160):finstr(txt160))
          write(txt12,fmt='(f12.4)') strMAD
          txt160 = 
     >    txt160(debstr(txt160):finstr(txt160))
     >    //'  &  '//txt12//'  \\'
        endif

        if(i.eq.1) then 
         if(typ2.eq.'intrinsic') then
          write(lunW,fmt='(10x,6('' & '',a12),'' & '',a,'' \\'')')  
     >    'Energy', 'Qz', ' e_z/pi', 'p_init', 'p_final', '|J_n|^2', 
     >    '\multicolumn{2}{c}{  |J_n|^2/ez  }'
          write(lunW,fmt='(10x,8('' & '',a12),'' \\'')')  
     >    '(GeV)', ' ', '(1e-6)', ' ',  ' ', '(1e-6)', 'zgoubi', 'MAD'
         else
c          write(lunW,fmt='(10x,6('' & '',a12),'' & '',a,'' \\'')')  
c     >    'Energy', 'Qz', ' zmax', 'p_init', 'p_final', '|J_n|^2', 
c     >    '\multicolumn{2}{c}{  |J_n|^2/ez  }'
c          write(lunW,fmt='(10x,8('' & '',a12),'' \\'')')  
c     >    '(GeV)', ' ', '  ', ' ',  ' ', '(1e-6)', 'zgoubi', 'MAD',
c     >    '(m)', '(rad)'
          write(lunW,fmt='(10x,5('' & '',a12),'' & '',a,'' \\'')')  
     >    'Energy', ' zmax', 'p_init', 'p_final', '|J_n|^2', 
     >    '\multicolumn{2}{c}{  |J_n|^2/ez  }'
          write(lunW,fmt='(10x,6('' & '',a12),'' \\'')')  
     >    '(GeV)', '  ', ' ',  ' ', '(1e-6)', 'zgoubi', 'MAD'
         endif
        endif
        write(lunW,fmt='(a)') txt160
      enddo

 123  continue

Create gnuplots 
      do i = 1, ifile
        cmmnd = 'cd '//
     >  drctry(i)(debstr(drctry(i)):finstr(drctry(i)))//' ; '//
     >  ' pwd ; '//
     >  dirTools(debstr(dirTools):finstr(dirTools))
     >  //'xing_dataTreatment/genePlots/genePlots'
        write(*,*) ' '
        write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
        call system(cmmnd)
      enddo

C Create ./Log/log.tex file
      cmmnd = 
     >dirTools(debstr(dirTools):finstr(dirTools))
     >//'xing_geneTexLog'//'/geneTexLog'
      write(*,*) ' '
      write(*,*) ' Pgm dataTreatment, execute : ',cmmnd
      call system(cmmnd)



      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'PGM DATATREATMENT COMPLETED ... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '
      stop

 99   continue
      write(*,*) 'Error open file scanSpinResonances.Out2'
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
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
      implicit double precision (a-h,o-z)
      LOGICAL STRCON
      CHARACTER STR*(*), STRIN*(*)
C     ------------------------------------------------------------------------
C     .TRUE. if the string STR contains the string STRIN with NCHAR characters
C     at least once.
C     IS = position of first occurence of STRIN in STR
C     ------------------------------------------------------------------------

      INTEGER DEBSTR,FINSTR

      II = 0
      DO 1 I = DEBSTR(STR), FINSTR(STR)
        II = II+1
        IF( STR(I:I+NCHAR-1) .EQ. STRIN ) THEN
          IS = II
          STRCON = .TRUE.
          RETURN
        ENDIF
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
