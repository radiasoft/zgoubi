C dataTreatment pgm treats data in M*Qz* directories. 
C This assumes prior execution of scanSpinResonances_launch.f in M*Qz*s' parent directory.
C dataTreatment is to be launched from the M*Qz*s' parent directory, it will 
C find (and get !) M*Qz*s' names from scanSpinResonances.Out2 storage file. 

      implicit double precision (a-h,o-z)

      character drctry(200)*132, zgDatFile(200)*132, txt132*132
      character typ2*12
      parameter(lunR=7,lunW=8)

      character cmmnd*300, fname*50
      integer debstr, finstr
      logical strcon

      parameter (pi = 4.*atan(1.),cl=2.99792458e8)

      character txtM*4, txtQz*15

      write(*,*) ' '
      write(*,*) '----------------------------------------------------'
      write(*,*) 'NOW RUNNING PGM COMPUTESTRENGTH... '
      write(*,*) '----------------------------------------------------'
      write(*,*) ' '

        fname = 
     >  './'
     >  //'zgoubi.res'
        open(unit=34,file=fname,err=991)
        read(34,fmt='(a)',end=102) txt132
        read(34,fmt='(a)',end=102) txt132
        read(34,*,end=102) boro

 22     continue
          read(34,fmt='(a)',end=102) txt132
          txt132 = txt132(debstr(txt132):132)
          if    (strcon(txt132,'''PARTICUL''',10,
     >                                           IS)) then
            read(34,*) am, q, G
          elseif(strcon(txt132,'''MULTIPOL''',10,
     >                                           IS)
     >    .and.  strcon(txt132,'SBEN',4,
     >                                  IS)
     >    .and.  strcon(txt132,'A1BF',4,
     >                                  IS)
     >                                            ) then
            read(34,fmt='(a)') txt132
            read(34,*) xl, ro, Bdip

          elseif(strcon(txt132,'''CAVITE''',8,
     >                                      IS)) then
            read(34,fmt='(a)') txt132
            read(34,*) circ, ah
            read(34,*) Vp, phis

          elseif(strcon(txt132,'DVCA02',6,
     >                                     IS)) then
C Get v-kick so to calculate zcoZgoubi = zcoMAD / vkick * (xl*b1/boro)
            read(34,fmt='(a)') txt132
            read(34,*) xl, ro, b1
C Normally 'SCALING' is used thus brho=1
            brho = 1000.d0
            vk = xl*b1/brho

          elseif(strcon(txt132,'***********************************',35,
     >                                  IS) ) then
            goto 102

          endif
        goto 22

 102    continue
        close(34)
     
C Get emittance and tune from tunesFromFai
        fname = 
     >  './'
     >  //'tunesFromFai.out'
        open(unit=34,file=fname)
        read(34,fmt='(a)',end=993) txt132
        read(34,fmt='(a)',end=993) txt132
        read(34,*,end=993) X0, Z0, XNU, ZNU, UXNU, UZNU, UX, UY, UZ
        call system(' pwd')
        write(*,*) i, X0, Z0, XNU, ZNU, UXNU, UZNU, UX, UY, UZ
        close(34)

C Get zco, vkick from MAD data
        fname = 
     >  './'
     >  //'scanSpinResonances.tmp'
        open(unit=34,file=fname)
        read(34,fmt='(a)',end=995) txt132
        read(34,*) 
     >  txtM, txtQz, aJn, aNn, dum, zco, vkick
        call system(' pwd')
        write(*,*) txtM, txtQz, aJn, aNn, dum, zco, vkick
        close(34)        

C Get Sz from avrgSzFromFai.out
        fname = 
     >  './'
     >  //'avrgSzFromFai.out'
        open(unit=34,file=fname)
        read(34,fmt='(a)',end=994) txt132
C First xing turns range
        read(34,fmt='(a)',end=994) txt132
        read(34,fmt='(a)',end=994) txt132
        read(34,fmt='(a)',end=994) txt132
        read(34,*,end=994) SzM,SzM2,Szmod,Szsig,SzImid
        read(34,fmt='(a)',end=994) txt132
C Last xing turns range
        read(34,fmt='(a)',end=994) txt132
        read(34,fmt='(a)',end=994) txt132
        read(34,fmt='(a)',end=994) txt132
        read(34,*,end=994) SzM,SzM2,Szmod,Szsig,SzFmid
        close(34)

C Get resonance style 
        fname = 
     >  '../'
     >  //'geneZGDat4Xing.data'
        open(unit=34,file=fname)
        read(34,fmt='(a)',end=996) typ2
        close(34)

        Bdip = Bdip / 10.d0  ! kG -> T
        am = am * 1d6        ! MeV -> eV
        p =  boro*cl/1.d3
        Etot = sqrt(p*p+am*am)
        DE = Vp * sin(phis)      
        alpha = G * DE / am / (2. * pi)
        rho = 1./Bdip
        R = circ / (2.*pi)
        Bdot = DE/(2.*pi*R*rho)
        ez = 2.d0*UY
        write(*,*) '  '
        write(*,*) ' Bdip(kG),circ,ah,Vp,phis,am,q,G,ez :'
     >          , Bdip,circ,ah,Vp,phis,am,q,G,ez
        write(*,*) '  '
        write(*,*) ' Present conditions : '
        write(*,*) '   xing speed alpha = dgamma/dtta = ',alpha
     >  , '    dB/dt  = ',Bdot


C Save results
        fname = 
     >  './'
     >  //'computeStrengths.out'
        open(unit=lunW,file=fname,err=992)
        pini = SzImid
        pfin = SzFmid
        A2 = -dlog((pfin/pini + 1.d0)/2.d0)
        if    (typ2(debstr(typ2):finstr(typ2)) .eq. 'intrinsic') then
          write(lunW,*) '% intrinsic ', 
     >    'E-tot, Qz, e_z/pi/1e-6, p_i, p_f, |J_n|^2/1e-6, |J_n|^2/ez'
          write(lunW,fmt='(1p,'' & '',g13.5,1x,'' & '',g13.5,1x,'' & '',
     >           g11.3,1x, 4('' & '',g12.4))') 
     >    Etot/1d9, uznu, ez*1.d6, pini, pfin, 
     >    A2*2.*alpha/pi*1.d6, A2*2.*alpha/pi/ez
          write(*,fmt='(1p,'' & '',g13.5,1x,'' & '',g13.5,1x,'' & '',
     >           g11.3,1x, 4('' & '',g12.4))') 
     >    Etot/1d9, uznu, ez*1.d6, pini, pfin, 
     >    A2*2.*alpha/pi*1.d6, A2*2.*alpha/pi/ez

        elseif(typ2(debstr(typ2):finstr(typ2)) .eq. 'imperfection') then
c          write(lunW,*) '% imperfection E-tot, Qz, zmax, p_i, p_f, |J_n|
c     >    ^2/1e-6, |J_n|^2/zmax^2 '
c          write(lunW,fmt='('' & '',g13.5,1x,'' & '',g13.5,1x,'' & '',
c     >           g11.3,1x, 4('' & '',f12.4))') 
c     >    Etot/1d9, uznu, zco*vk/vkick, pini, pfin, A2*2.*alpha/pi*1.d6
c     >    ,A2*2.*alpha/pi/(zco*vk/vkick)**2
          write(lunW,*) '% imperfection E-tot, zmax, p_i, p_f, |J_n|
     >    ^2/1e-6, |J_n|^2/zmax^2 '
          write(lunW,fmt='('' & '',g13.5,1x,'' & '',
     >           g11.3,1x, 4('' & '',f12.4))') 
     >    Etot/1d9, zco*vk/vkick, pini, pfin, A2*2.*alpha/pi*1.d6
     >    ,A2*2.*alpha/pi/(zco*vk/vkick)**2
        endif
      stop

 99   stop ' Pgm computeStrengths : error open scanSpinResonances.Out2'
 991  stop ' Pgm computeStrengths : error open zgoubi.res'
 992  stop ' Pgm computeStrengths : error open computeStrengths.out'
 993  stop ' Pgm computeStrengths : no data in tunesFromFai.out !'
 994  stop ' Pgm computeStrengths : no data in avrgSzFromFai.out !'
 995  stop ' Pgm computeStrengths : no data in scanSpinResonances.tmp'
 996  stop ' Pgm computeStrengths : no data in geneZGDat4Xing.data'
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


