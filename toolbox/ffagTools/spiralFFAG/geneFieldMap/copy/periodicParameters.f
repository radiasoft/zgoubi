C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C----- PLOT SPECTRUM     
      implicit double precision (a-h,o-z)
      CHARACTER txt132*132, txt9*9
      dimension xco(99), xpco(99), zco(99), zpco(99), dpco(99), tof(99)
      dimension path(99), alpha(99), eta(99)
      dimension betx(99), alfx(99), Dx(99)
      dimension betz(99), alfz(99), Dz(99)
      dimension xnu(99), znu(99)
      LOGICAL IDLUNI
      LOGICAL STRCON,  first
      CHARACTER LET*1
      data first / .true. /
      data lunRKx / 9 /

c---------------- Hypothesis data : 
      call readat(lunIn,nCell)
      write(6,*) ' Input data read from lattice.data file, ok !'
c---------------- END HYPOTHESIS 

C Open zgoubi.res
      IF (IDLUNI(lunR)) THEN
        OPEN(UNIT=lunR,FILE='zgoubi.res',ERR=799)
      ELSE
        GOTO 798
      ENDIF

C Search zgoubi.res for co coordinates
        read(lunR,fmt='(a)',end=62) txt132
        read(lunR,fmt='(a)',end=62) txt132
        read(lunR,fmt='(a)',end=62) txt132
        read(lunR,fmt=*,end=62) objk
        read(lunR,fmt=*,end=62) nco
        write(*,*) ' '
        write(*,*) '----------------------------'
        write(*,*) 'Computation of periodic parameters, based on ',
     >    'the following ',nco,' closed orbits :'
        do ico = 1, nco
          read(lunR,*,end=62)
     >    xco(ico),xpco(ico),zco(ico),zpco(ico),dum,dpco(ico)
          write(*,fmt='(1p,5g14.6,i6)')
     >    xco(ico),xpco(ico),zco(ico),zpco(ico),dpco(ico),ico
        enddo
        write(*,*) ' '

C Search zgoubi.res for "'MARKER' #END"
      first = .true.
 60   continue
        read(lunR,fmt='(a)',end=62) txt132
        write(lunW,*) txt132   
        if(.not.strcon(txt132,'#END',4,
     >                                 IS)) then
          goto 60
        else
          if(first) then 
             first = .false. 
             goto 60
          endif
        endif

        jjco = 1
 61     continue

        read(lunR,fmt='(a)',end=62) txt132
c        write(lunW,*) txt132   

        if(strcon(txt132,'flight',6,
     >                              IS)) then
            backspace(lunR)
            backspace(lunR)
            read(lunR,fmt='(a)',end=62) txt132
c               write(*,*) txt132
            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1,F2,F3,F4,F5,F6,  I
 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F11.3,1X,I6)
            path(jjco) = f6
c  read tof of traj 1 
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tof(jjco)
c            write(*,*) ' Path, TOF of traj. ',I,' : ',path(jjco),' cm',
c     >             tof(jjco),' mu_s'
c  read path length of traj 2 and 3 
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1p,F2,F3,F4,F5,F6p,  I
c            write(*,*) ' Path length of traj. # ',I,' is ',F6p,' cm'
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tofp
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1m,F2,F3,F4,F5,F6m,  I
c            write(*,*) ' Path length of traj. # ',I,' is ',F6m,' cm'
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tofm
            
            dpp = (f1p-f1m) / f1
            alpha(jjco) = (f6p-f6m)/f6 / dpp
            eta(jjco) = (tofp-tofm)/tof(jjco) / dpp

c            write(*,*)'jjco=',jjco,dpp,f1p,f1m,f1,alpha(jjco),eta(jjco)
c            write(*,*)'jjco=',jjco,dpp,f6p,f6m,f6,alpha(jjco),eta(jjco)

            jjco = jjco + 1
        endif

      goto 61
 62   close(lunR)

c      call system(
c     > 'cat periodicParameters.out >> periodicParameters.out_cat')
            IF (IDLUNI(LUNW)) THEN
              OPEN(UNIT=LUNW,FILE='periodicParameters.out',ERR=699)
            ELSE
              GOTO 698
            ENDIF
c      WRITE(LUNW,*) ' % path(cm),tof(mu_s),alpha,eta (ncell=',ncell,')',
c     >   '   K=',xK,'  xi=',xiDeg,' deg.'

          open(unit=34,file='tempKXi.dum')
          read(34,*) xK,xiDeg
          close(34)

      write(*,*)
      do ico = 1, jjco-1
        WRITE(*,179) path(ico)*nCell,tof(ico)*nCell,
     >     alpha(ico)/10,eta(ico)/10,ico
        WRITE(LUNW,179) path(ico)*nCell,tof(ico)*nCell,
     >     alpha(ico)/10,eta(ico)/10,ico, xK,xiDeg,
     >    ' % path(cm),tof(mu_s),alpha,eta (ncell=',ncell,
     >    '),  K,  xi (deg.)'
c 179    FORMAT(1P,e14.6,e16.8,2f10.4,i6,)
 179    FORMAT(1P,e14.6,e16.8,2f10.4,i6,2e12.4,a,i6,a)
      enddo

      write(6,*) ' Job ended.  Went on well it seems...'
      goto 99

 698  WRITE(6,*) '*** Problem : No idle unit for periodicParameters.out'
      GOTO 99
 699  WRITE(6,*) '*** Problem at OPEN periodicParameters.out '
      GOTO 99

 798  WRITE(6,*) '*** Problem : No idle unit for open zgoubi.res '
      GOTO 99
 799  WRITE(6,*) '*** Problem at OPEN zgoubi.res '
      GOTO 99

 99   continue
      stop
      end
      FUNCTION IDLUNI(LN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI

      LOGICAL OPN

      I = 20
 1    CONTINUE
        INQUIRE(UNIT=I,ERR=99,IOSTAT=IOS,OPENED=OPN)
        I = I+1
        IF(I .EQ. 100) GOTO 99
        IF(OPN) GOTO 1
        IF(IOS .GT. 0) GOTO 1
      
      LN = I-1
      IDLUNI = .TRUE.
      RETURN

 99   CONTINUE
      LN = 0
      IDLUNI = .FALSE.
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

      subroutine readat(lunIn,
     >                        nCell)
      implicit double precision (a-h, o-z)
      logical strcon, ok
      character*132 txt132
      integer debstr

      open(unit=lunIn,file='lattice.data')

      read(lunIn,fmt='(a)',end=99,err=98)  txt132  !! reads comment line. Next line is to be K
      write(*,*) txt132
 
 1    continue

      read(lunIn,fmt='(a)',end=99,err=98)  txt132
      ideb = debstr(txt132)
      txt132 = txt132(ideb:132)
c      write(*,*) txt132

      if (STRCON(txt132(1:6),'K',1,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) aK
        write(*,*) ' aK = ', aK
      elseif (STRCON(txt132(1:6),'xiDeg',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) xi
        write(*,*) ' xiDeg = ', xi
      elseif (STRCON(txt132(1:6),'nCell',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) nCell
        write(*,*) ' nCell = ', nCell
      elseif (STRCON(txt132(1:6),'r0',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) r2
        write(*,*) ' r2 = ', r2
      elseif (STRCON(txt132(1:6),'pf',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) pf
        write(*,*) ' pf = ', pf
      elseif (STRCON(txt132(1:6),'gap',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) gap
      elseif (STRCON(txt132(1:6),'kappa',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) akappa
        write(*,*) ' kappa = ', akappa
      elseif (STRCON(txt132(1:6),'T1',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) T1
      elseif (STRCON(txt132(1:6),'T2',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) T2
      elseif (STRCON(txt132(1:6),'nCO',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) nCO
        write(*,*) ' nCO = ', nCO
      elseif (STRCON(txt132(1:6),'map',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) map
      endif

      goto 1

 99   continue
      close(lunIn)
      write(*,*) ' Read lattice.data ended upon eof'
      return
 98   continue
      close(lunIn)
      write(*,*) ' Read lattice.data ended upon read error'
      return
      end

      subroutine gotoEnd(lun,
     >                       nLine)
      implicit double precision (a-h,o-z)
      character*132 txt132

      write(6,*) '  Now go to End of file lun=',lun

      nLine = 0
 10   read(lun,*,end=98,err=99) txt132
        nLine = nLine + 1 
        write(*,*) ' gotoEnd, ', nLine, txt132
        goto 10

      return

 98   write(6,*) '  * Now at end of file !! *'
      return
 99   write(6,*) '  * gotoEnd terminated upon read error *'
      return
      end
