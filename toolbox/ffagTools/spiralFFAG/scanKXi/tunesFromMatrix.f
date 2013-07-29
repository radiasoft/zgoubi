C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C----- PLOT SPECTRUM     
      implicit double precision (a-h,o-z)
      CHARACTER txt132*132, txt9*9
      dimension xco(99), xpco(99), zco(99), zpco(99), dpco(99), tof(99)
      dimension path(99)
      dimension betx(99), alfx(99), Dx(99)
      dimension betz(99), alfz(99), Dz(99)
      dimension xnu(99), znu(99)
      LOGICAL IDLUNI
      LOGICAL STRCON, EXS
      CHARACTER LET*1
      CHARACTER*14 txtxnu, txtznu

C Just a temp storage for passing K, Xi :
      INQUIRE(FILE='tempKXi.dum',exist=EXS)
      if(exs) then
        open(unit=34,file='tempKXi.dum')
        read(34,*,err=10,end=10)  xK, xiDeg
        close(34)
      else
        xK= 0.
        xiDeg = 0. 
      endif
 10   continue

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
        kobj = objk
        kobj2 = nint(100.*objk) - 100*kobj
        if(kobj2.eq.0) then
          nCO = 1
        else
          nCO = kobj2
        endif
        read(lunR,fmt='(a)',end=62) txt132
        do ico = 1, nco
          read(lunR,*,end=62)
     >    xco(ico),xpco(ico),zco(ico),zpco(ico),dum,dpco(ico)
        enddo

C Search zgoubi.res for MATRIX outputs
      iico = 1
      jjco = 1
 61   continue

        read(lunR,fmt='(a)',end=62) txt132
c        write(lunW,*) txt132   

        if(strcon(txt132,'Beam  matrix  (beta',19,
     >                                            IS)) then 
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,*,end=62) betx(iico),dum,dum,dum,dum,Dx(iico)
          read(lunR,*,end=62) alfx(iico),dum,dum,dum,dum,Dpx
          alfx(iico) = -alfx(iico)
          read(lunR,*,end=62) dum,dum,betz(iico)
          read(lunR,*,end=62) dum,dum,alfz(iico)
          alfz(iico) = -alfz(iico)
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,*,end=62) 
     >          txt9,txt9, txtxnu, txt9,txt9, txtznu
          if(strcon(txtxnu,'undefined',9,
     >                                   IS)) then 
            xnu(iico) = -1.d0
          else
            read(txtxnu,*)  xnu(iico)
          endif
          if(strcon(txtznu,'undefined',9,
     >                                   IS)) then 
            znu(iico) = -1.d0            
          else
            read(txtznu,*)  znu(iico)
          endif

          iico = iico + 1
        endif

        if(iico .gt. 1) then      !! means MATRIX has already been seen
          if(strcon(txt132,'Time of flight',14,
     >                                         IS)) then 

c  read path length of traj 1 (modulo 11)
            backspace(lunR)
            backspace(lunR)
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1,F2,F3,F4,F5,F6,  I
 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F11.3,1X,I6)
            path(jjco) = f6
c  read tof of traj 1 (modulo 11)
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tof(jjco)
            write(*,*) ' Path, TOF of traj. ',I,' : ',path(jjco),' cm',
     >             tof(jjco),' mu_s'

c  jump to traj. 10 (modulo 11)
            do jj = 1, 16
              read(lunR,fmt='(a)',end=62) txt132
            enddo

c  read path length of traj 10 and 11 (modulo 11)
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

            jjco = jjco + 1
          endif
        endif

      goto 61
 62   close(lunR)

      call system('cat tunesFromMatrix.out >> tunesFromMatrix.out_old ; 
     >                   rm -f tunesFromMatrix.out ')
            IF (IDLUNI(LUNW)) THEN
              OPEN(UNIT=LUNW,FILE='tunesFromMatrix.out',ERR=699)
            ELSE
              GOTO 698
            ENDIF
      WRITE(LUNW,*) 
     >' xco,xpcoXNU,ZNU,alfx,betx,alfz,betz,Dx,p/p0,tof(s),kt,npts,K,xi'

      do ico = 1, nco
        WRITE(LUNW,179) xco(ico),xpco(ico), 
     >  XNU(ico),ZNU(ico),alfx(ico),betx(ico),alfz(ico),betz(ico),
     >  Dx(ico),dpco(ico), tof(ico), KT,NPTS, xK, xiDeg
 179    FORMAT(1P,10G14.6,G18.10,2I4,2G12.4)
      enddo

      write(6,*) ' Job ended.  Went on well it seems...'
      goto 99

 698  WRITE(6,*) ' *** Problem : No idle unit for tunesFromMatrix.out '
      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN tunesFromMatrix.out '
      GOTO 99

 798  WRITE(6,*) ' *** Problem : No idle unit for open zgoubi.res '
      GOTO 99
 799  WRITE(6,*) ' *** Problem at OPEN zgoubi.res '
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
