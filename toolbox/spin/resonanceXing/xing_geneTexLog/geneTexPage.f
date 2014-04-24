      character txt132*132, txtfig*20
      parameter (lunR=7,lunW=8,lunMQ=9)
      character txtM*4, txtQ*14
      integer debstr, finstr
      
      open(unit=lunR,file='geneTexLog.out')
      read(lunR,*) txtfig
      close(lunR)

      txtfig = txtfig(debstr(txtfig):finstr(txtfig))
      i = 1
 1    continue
        if(txtfig(i:i) .eq. 'Q') then
          txtQ = txtfig(i+2:finstr(txtfig))
          if(txtQ .eq. '0') txtQ = ' '
          txtM = txtfig(debstr(txtfig)+1:i-1)
        else
          i = i + 1
          goto 1
        endif

      open(unit=lunW,file='logPage.tex')

      txt132 = '\clearpage'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = '\begin{figure}[h]'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = '\addcontentsline{toc}{subsubsection}{\hspace{10ex} '
     >  //' $\gamma G = '//txtM//' '//txtQ//'$ }'
      write(lunW,*) txt132
      txt132 = '{ $\gamma G = '//txtM//' '//txtQ//'$ }'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = '\begin{center}'
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=14cm,height=5cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xingFull.eps}'
      write(lunW,*) txt132
      txt132 = '  \caption{ \label{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))//'_xingFull}'
      write(lunW,*) txt132
      txt132 = '$S_y$ versus kinetic energy. '
      write(lunW,*) txt132
      txt132 = ' }'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132

      txt132 = '\hspace{-.5cm}'
      write(lunW,*) txt132
      txt132 = '\mbox{'
      write(lunW,*) txt132
      txt132 = '\begin{minipage}[c]{0.49\linewidth}'
      write(lunW,*) txt132
      txt132 = '\centering'
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=8cm,height=4cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xingInit.eps}'
      write(lunW,*) txt132
      txt132 = '  \caption{ \label{' 
     >  //txtfig(debstr(txtfig):finstr(txtfig))//'_xingInit}'
      write(lunW,*) txt132
      txt132 = 'Zoom on initial $S_y$. '
      write(lunW,*) txt132
      txt132 = ' }'
      write(lunW,*) txt132
      txt132 = '\end{minipage}\hspace{0.05\linewidth}'
      write(lunW,*) txt132
      txt132 = '\begin{minipage}[c]{0.49\linewidth}'
      write(lunW,*) txt132
      txt132 = '  \centering'
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=8cm,height=4cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xingFinal.eps}'
      write(lunW,*) txt132
      txt132 = '  \caption{ \label{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))//'_xingFinal}'
      write(lunW,*) txt132
      txt132 = 'Zoom on final $S_y$. '
      write(lunW,*) txt132
      txt132 = ' }'
      write(lunW,*) txt132
      txt132 = '\end{minipage}\hspace{0.05\linewidth}'
      write(lunW,*) txt132
      txt132 = '}'   !mbox
      write(lunW,*) txt132

      txt132 = ' '
      write(lunW,*) txt132
      txt132 = ' \end{center}'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = ' \begin{center}'
      write(lunW,*) txt132

      txt132 = '\mbox{'
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=8cm,height=4cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xing-xxp.eps}'
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=8cm,height=4cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xing-yyp.eps}'
      write(lunW,*) txt132
      txt132 = '}'   !mbox
      write(lunW,*) txt132

      txt132 = '  \caption{ \label{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))//'_xing-xxp}'
      write(lunW,*) txt132
      txt132 = 'x-x'' and z-z''. '
      write(lunW,*) txt132
      txt132 = ' }'
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = ' '
      write(lunW,*) txt132
      txt132 = '\includegraphics*[width=8cm,height=4cm]{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))
     >  //'_xing-ldp.eps}'
      write(lunW,*) txt132
      txt132 = '  \caption{ \label{'
     >  //txtfig(debstr(txtfig):finstr(txtfig))//'_xing-dp}'
      write(lunW,*) txt132
      txt132 = 'dp-phase. '
      write(lunW,*) txt132
      txt132 = ' }'
      write(lunW,*) txt132
      txt132 = '  \end{center}'
      write(lunW,*) txt132
      txt132 = '\end{figure}'
      write(lunW,*) txt132

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

      
