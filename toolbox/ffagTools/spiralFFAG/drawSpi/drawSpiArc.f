C      implicit double precision (a-h, o-z)
      data lunIn, lunOut, lunres / 7, 8, 10 /
      data zero / 0. /
      pi = 4.d0 * atan(1.d0)
      deg2rd = pi / 180.d0
      c = 2.99792458d8
      am = 938.27231d6


      r2 = 4.1
      r1 = 2.2
      xi = 47.1
 
      write(6,fmt='(a,F10.6,a,2f10.6,a)') 
     > ' xi,  r1, r2 : ', xi,' (deg.)  ',  r1, r2,' (m)'

      r0 = (r2+r1)/2.d0
      dr = (r2-r1)/2.d0

      write(lunRes,*)  ' '
      write(lunRes,*) 'It also results from that all : '
      write(lunRes,*) ' drift_inj = ',2.d0*pi*(r0-dr)*(1.d0-pf)/nCell,
     >' m,    drift_xtr = ',2.d0*pi*(r0+dr)*(1.d0-pf)/nCell,' m'  

      close(lunRes)

      b = 1.d0 / tan(xi*pi/180.d0)

        tta0 = 0.
C Limit angles corresponding to limit radii r1, r2
        tta1 =  log(r1/r0) / b
        tta2 =  log(r2/r0) / b

          write(*,*) '  ttaMin, max : ', tta1, tta2
c         pause

Compute first spiral EFB
        do tta = tta1, tta2*1.0001, (tta2-tta1)/20.d0
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0)  
           y = r *sin(tta + tta0)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(*,fmt='(2F10.6,F6.1)') r, tta
           write(88,fmt='(2F10.6,F6.1)') x, y, zero
        enddo
c      enddo
      
      stop
      end
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
C      implicit double precision (a-h, o-z)
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
C      implicit double precision (a-h, o-z)
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
C      implicit double precision (a-h, o-z)
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

      subroutine readat(lunIn,nCell, aK, xi, pf, r2, gap, T1, T2)
C      implicit double precision (a-h, o-z)
      logical strcon
      character*132 txt132
      integer debstr

      open(unit=lunIn,file='drawSpi.data')

      read(lunIn,fmt='(a)',end=99)  txt132  !! reads comment line. Next line is to be K

      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) aK
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) xi
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) nCell
      endif 
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) r2
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) pf
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) gap
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) T1
      endif
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) T2
      endif

 99   continue
      close(lunIn)
      return
      end














