      character txt132*132, txt80*80, txtfig*20, txtM*4, txtQ*20
      parameter (lunW=8,lunR=7)
      integer debstr, finstr

      open(unit=lunW,file='./Log_Chroma/log_Chroma.tex',status='old')

      call gotoEnd(lunW,
     >                  line)

      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\clearpage'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
c      txt80='\addtocontents{toc}{\protect\vspace{6pt}}%'
c      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\addcontentsline{toc}{section}{Appendix}%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\appendix'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\section{Zgoubi data file (top of)}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='{\tiny'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\begin{verbatim}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))

      open(unit=lunR,file='zgoubi.dat',status='old')

      do i = 1, 100
        read(lunR,fmt='(a)',err=11,end=11) txt132
        write(lunW,*) txt132
      enddo
 11   continue
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\end{verbatim}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='} %\footnotesize'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))

      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\end{document}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\smallskip'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))

      close(lunR)
      close(lunW)
      stop
      end
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
      write(*,*) 'Pgm geneTexEnd.  Error reading in lunW ' 
      write(*,*) 'line = ',line
      write(*,*) '************* '
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

