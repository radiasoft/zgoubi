      character*132 txt132
      character*100 cmmnd
      parameter(lunR=7, lunW=8)
      integer debstr, finstr
      data txt132 / '\smallskip'  /

      open(unit=lunR,file='gnuplot_alfa.Out',status='old')
      read(lunR,*,end=10,err=10) txt132
      read(lunR,*,end=10,err=10) aL0, gtr, alfa, alfa1, alfa2
 10   continue
      close(lunR)
      open(unit=lunR,file='gnuplot_eta.Out',status='old')
      read(lunR,*,end=11,err=11) txt132
      read(lunR,*,end=11,err=11) T0, eta, eta1, eta2 
 11   continue
      close(lunR)      
      open(unit=lunR,file='gnuplot_dQdp.Out',status='old')
      read(lunR,*,end=12,err=12) Qx0,Qx1,Qx2 !,Qx3
      read(lunR,*,end=12,err=12) Qy0,Qy1,Qy2 !,Qy3
 12   continue
      close(lunR)

C Open ./Log_Chroma/log_Chroma.tex and write body
      open(unit=lunW,file='./Log_Chroma/log_Chroma.tex',status='old')
      call gotoEnd(lunW,
     >                  line)

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132='\section{Orbits, dispersion}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\begin{center}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\includegraphics*['
     >//'width=17cm,height=7cm]{'
     >  //'gnuplot_CO.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigCOxy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Closed orbits, H and V. '
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
      txt132 = '\includegraphics*['
     >//'width=17cm,height=7cm]{'
     >  //'gnuplot_DxDy.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigDxDy}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Dispersion, H and V (from orbit '
     >//'difference - getDiffFromFai).'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' \end{center}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      txt132 = '\end{figure}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))


      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
C      txt132 = ' \clearpage'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

C      txt132='\section{Orbits, dispersion, $\alpha$, $\eta$, Chromaticity}'
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132='\section{Optical functions}'
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
      txt132 = 'Optical functions, H and V '
     >//'(from $\beta$-matrix transport - betaFromMatrix).'
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
      txt132 = 'Dispersion, H and V '
     >//'(from $\beta$-matrix transport - betaFromMatrix).'
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
      txt132 = 'Phase advance, H and V '
     >//'(from $\beta$-matrix transport - betaFromMatrix).'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\end{figure}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = '\clearpage'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132='\section{$\alpha$, $\eta$, tunes, chromaticities}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\begin{center}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\includegraphics*['
     >//'width=14cm,height=5cm]{'
     >  //'gnuplot_alfa.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigAlfa}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Momentum dependence of orbit length.'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\vspace{2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) '$L_0 = ',aL0
     >,'$, $\gamma_{tr} = ',gtr
     >,'$, $\alpha_0 = ',alfa,'$,  $\alpha_1 = ',alfa1
     >,'$, $\alpha_2 = ',alfa2,'$'
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' \vspace{4ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\rule{100mm}{0.2mm}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' \end{center}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      txt132 = '\end{figure}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' ~ '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' ~ '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))


c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\begin{figure}[h]'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\begin{center}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\includegraphics*['
     >//'width=14cm,height=5cm]{'
     >  //'gnuplot_eta.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigEta}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Momentum dependence of revolution time.'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\vspace{2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) '$T_0 = ',T0
     >,'$, $\eta = ',eta,'$,  $\eta_1 = ',eta1 
     >,'$, $\eta_2 = ',eta2,'$'
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' \vspace{4ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\rule{100mm}{0.2mm}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' \end{center}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      txt132 = '\end{figure}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))




c     txt132 = '\clearpage'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\begin{figure}[h]'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\begin{center}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\includegraphics*['
     >//'width=14cm,height=5cm]{'
     >  //'gnuplot_Chroma.eps}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\vspace{-2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '  \caption{ \label{'
     >//'FigdQdp}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = 'Momentum dependence of tunes, and '
     >//' matching polynomials $Q(\delta) = Q_0 + Q''\delta '
     >//' + Q'''' \delta^2 + Q'''''' \delta^3$. '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' }'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '\vspace{2ex}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      write(lunW,*) 'Horizontal~:  $Q_x = ',Qx0,'$, $Q_x'' = ',Qx1
     >,'$, $Q_x'''' = ',Qx2,'$' 
C     >,'$, $Q_x'''' = ',Qx2,'$,  $Q_x'''''' = ',Qx3,'$' 
      txt132 = '    '
      write(lunW,*) txt132(1:3)
      write(lunW,*) 'Vertical~:  $Q_y = ',Qy0,'$, $Q_y'' = ',Qy1
     >,'$, $Q_y'''' = ',Qy2,'$' 
C     >,'$, $Q_y'''' = ',Qy2,'$,  $Q_y'''''' = ',Qy3,'$' 
      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' \vspace{5ex}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\rule{100mm}{0.2mm}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = 'MONITORING'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' \vspace{5ex}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = '\includegraphics*['
c     >//'width=14cm,height=5cm]{'
c     >//'gnuplot_xxp.eps}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\vspace{-2ex}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '  \caption{ \label{'
c     >  //'Figxxp}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' Horizontal phase-space motions used for computing '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' Qx values as shown in Fig.~\ref{FigdQdp}. ' 
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' }'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = '\vspace{6ex}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\includegraphics*['
c     >//'width=14cm,height=5cm]{'
c     >//'gnuplot_yyp.eps}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '\vspace{-2ex}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = '  \caption{ \label{'
c     >  //'Figyyp}'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' Vertical phase-space motions used for computing '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' Qy values as shown in Fig.~\ref{FigdQdp}. ' 
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' }'
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
c      txt132 = ' '
c      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = ' '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = ' \end{center}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\end{figure}'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\smallskip'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))
      txt132 = '\smallskip'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      close(lunW)
      
cC End section of geneTexlog.tex
c      cmmnd = 
c     >'~/zgoubi/SVN/current/toolbox/spin/xing_geneTexLog/geneTexEnd'
c      call system(cmmnd)

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
      if(length .ne. 0) then
1       CONTINUE
          DEBSTR=DEBSTR+1
C          IF(DEBSTR .EQ. LENGTH) RETURN
C          IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') GOTO 1
          IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
            IF(DEBSTR .EQ. LENGTH) THEN
              DEBSTR = 1
              RETURN
            ELSE
              GOTO 1
            ENDIF
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
        IF(FINSTR .EQ. 1) THEN
           RETURN
        ENDIF
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END

