      character txt80*80, txtfig*20, txtM*4, txtQ*20
      parameter (lunW=8)
      integer debstr, finstr

      open(unit=lunW,file='./Log_Chroma/log_Chroma.tex')

      txt80='\documentclass[11pt]{article}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{draftcopy}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%\usepackage[draft]{graphicx}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{graphicx}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{wrapfig}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{amssymb}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{lscape}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\usepackage{times}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\oddsidemargin =-0.1in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\evensidemargin=-0.1in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\textwidth=6.8in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\textheight=10.2in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\topmargin=-.8in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\footskip=0in'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\Br}{\ensuremath{B\!\rho}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\bull}{\ensuremath{\bullet~}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\nib}{\noindent \ensuremath{\bullet~}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\com}{{center of mass}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\hbrk}{\hfill \break}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\lab}{{laboratory frame}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\MC}{Monte~Carlo}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\rms}{\ensuremath{rms}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\wrt}{{with respect to}}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\rfrncA}{\rm  }'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\rfrncB}{\rm }'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\newcommand{\rfrncC}{\rm }'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\pagestyle{myheadings}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\markboth{\small \rfrncA \rfrncB  \rfrncC \hfill}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='         {\small \rfrncA \rfrncB  \rfrncC \hfill}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\pagestyle{headings}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\begin{document}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\thispagestyle{myheadings}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='~ ~ ~ '
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80=' ~ ~ ~ \hfill F. Meot '
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\begin{center}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\bf \large '
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='Output from "SearchChroma" data treatment'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\end{center}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
C      txt80='\vfill'
C      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\tableofcontents'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\bigskip '
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\rule{100mm}{0.1mm}'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='%'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))
      txt80='\smallskip'
      write(lunW,*) txt80(debstr(txt80):finstr(txt80))

C      txt132 = '\clearpage'
C      write(lunW,*) txt132

      close(lunW)
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
