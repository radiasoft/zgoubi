C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C----- PLOT SPECTRUM     
      implicit double precision (a-h,o-z)
      CHARACTER txt132*132, txt9*9
      dimension xco(4999),xpco(4999),zco(4999),zpco(4999)
      dimension dpco(4999),tof(4999)
      dimension path(4999)
      dimension betx(4999), alfx(4999), Dx(4999), Dpx(4999)
      dimension betz(4999), alfz(4999), Dz(4999), Dpz(4999)
      dimension phix(4999), phiz(4999), numlm(4999)
      LOGICAL IDLUNI
      LOGICAL STRCON, EXS
      CHARACTER LET*1
      logical okpart
      CHARACTER KEYL12*50, txt1*1

      INTEGER DEBSTR,FINSTR
      data okpart / .false. /

      logical part
C EMMA
      data am, qe /  0.51099892, 1.602176487e-19  /

C Open getScumFromRes.in
      IF (IDLUNI(lunIn)) THEN
        OPEN(UNIT=lunIn,FILE='getScumFromRes.in',ERR=797)
C        read(lunIn,fmt='(a)') KEYL12
        read(lunIn,*) KEYL12
C        write(*,*) ' Kley and labels to be found : ',KEYL12,
C     >' ---',KEYL12(debstr(KEYL12):finstr(KEYL12)),'---'
      ELSE
        GOTO 797
      ENDIF
      close(lunIn)

C Open getScumFromRes.out
      IF (IDLUNI(lunW)) THEN
        OPEN(UNIT=lunW,FILE='getScumFromRes.out',ERR=796)
      ELSE
        GOTO 796
      ENDIF

C Open zgoubi.res
      IF (IDLUNI(lunR)) THEN
        OPEN(UNIT=lunR,FILE='zgoubi.res',ERR=799)
      ELSE
        GOTO 798
      ENDIF

C Search zgoubi.res for KEYL12
      lng = len(KEYL12(debstr(KEYL12):finstr(KEYL12)))
C        write(*,*) 'Length of string : ', lng

      noc = 1
 61   continue

        read(lunR,fmt='(a132)',end=62) txt132
C        write(*,fmt='(a132)') txt132

        if(strcon(txt132,KEYL12(debstr(KEYL12):finstr(KEYL12)),lng,
     >                             IS)) then 

C        write(*,fmt='(a132)') txt132

 611      continue
          read(lunR,fmt='(a132)',end=62) txt132

          if(strcon(txt132,'Cumulative length of opti',25,
     >                                            IS)) then 

            read(txt132,*) txt1, txt1, txt1, txt1, txt1, txt1, dist
C            write(*,fmt='(i6,1p,2x,e16.8,2x,a)') noc, dist,
            write(lunW,fmt='(i6,2x,f12.6,2x,a)') noc, dist,
     >      KEYL12(debstr(KEYL12):finstr(KEYL12))
            noc = noc + 1

          elseif(strcon(txt132,'***********************************',35,
     >                                            IS)) then 
            goto 61

          else
            goto 611
          endif
        endif

      goto 61
 62   close(lunR)

      write(6,*) ' Job ended.  Went on well it seems...'
      goto 99

 698  WRITE(6,*) ' *** Problem : No idle unit for getScumFromRes.out '
      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN getScumFromRes.out '
      GOTO 99

 796  WRITE(6,*) ' *** Problem opening getScumFromRes.out'
      GOTO 99

 797  WRITE(6,*) ' *** Problem opening getScumFromRes.in '
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
