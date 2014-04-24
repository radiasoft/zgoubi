      character*132 txt132
      character*100 cmmnd
      parameter(lunR=7, lunW=8)
      INTEGER DEBSTR,FINSTR

C Open ./Log_CO/log_CO.tex and write body
      open(unit=lunW,file='./Log_CO/log_CO.tex',status='old')
      call gotoEnd(lunW,
     >                  line)
      
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\centering'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\includegraphics*[width=17cm,height=7cm]{'
     >  //'gnuplot_CO.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigCOxy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Closed orbits, H and V, from zgoubi.fai. '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\includegraphics*[width=17cm,height=7cm]{'
     >  //'gnuplot_DxDy.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigDxDy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Dispersion, H and V (getDiffFromFai).'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\end{figure}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' \clearpage'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\centering'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\includegraphics*[width=17cm,height=6cm]{'
     >  //'gnuplot_betxy.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'Figbetxy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Optical functions, H and V (betaFromMatrix).'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\includegraphics*[width=17cm,height=6cm]{'
     >  //'gnuplot_Dxy.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigDxy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Dispersion, H and V (betaFromMatrix). '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      txt132 = ' ~'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\includegraphics*[width=17cm,height=6cm]{'
     >  //'gnuplot_muxy.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'Figmuxy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Phase advance, H and V (betaFromMatrix). '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\end{figure}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132=' '
      write(lunW,*) txt132
      txt132=' '
      write(lunW,*) txt132



      close(lunW)
      stop
      end

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
      FUNCTION DEBSTR(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER * (*) STR
C     -----------------------------------
C     Renvoie dans DEBSTR le rang du
C     premier caractere non-blanc de STR.
C     -----------------------------------
      DEBSTR=0
      LENGTH=LEN(STR)+1
1     CONTINUE
         DEBSTR=DEBSTR+1
         IF(DEBSTR.GE. LENGTH) RETURN
         IF (STR(DEBSTR:DEBSTR).EQ. ' ') GOTO 1
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
      subroutine gotoEnd(lunW,
     >                        line)
      character*132 txt132
      line = 0
 1    continue
        read(lunW,*,err=99,end=10) txt132
        line = line + 1
        goto 1

 10   continue
      backspace(lunW)
      return

 99   continue
      write(*,*) ' '
      write(*,*) '************* '
      write(*,*) 'Pgm geneTexBody.  Error reading in lunW ' 
      write(*,*) 'line = ',line
      write(*,*) '************* '
      write(*,*) ' '
      stop 
      end

