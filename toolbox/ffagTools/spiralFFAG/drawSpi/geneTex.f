      data lunIn, lunOut / 11, 12 /
      character*132 txt132, prtcl*6
      integer debstr, finstr
      data prtcl / 'proton' /

      pi = 4. * atan(1.)
      r2d = 180./pi

      write(6,*) 
      write(6,*) '--------------------------------------------'
      write(6,*) ' geneTex procedure now going on...' 

      open(unit=lunOut,file='log.tex')

      write(lunOut,*) '\\documentclass[11pt]{article}'
      write(lunOut,*) '\\usepackage{graphicx}'
      write(lunOut,*) '\\usepackage{wrapfig}'
      write(lunOut,*) '\\usepackage{amssymb}'
      write(lunOut,*) '\\usepackage{lscape}'
      write(lunOut,*) '\\usepackage{times}'
      write(lunOut,*) '\\usepackage{color}'
      write(lunOut,*) ''
      write(lunOut,*) '\\oddsidemargin=0.15in'
      write(lunOut,*) '\\evensidemargin=0.15in'
      write(lunOut,*) '\\textwidth=6.6in'
      write(lunOut,*) '\\textheight=9.3in'
      write(lunOut,*) '\\topmargin=0.in'
      write(lunOut,*) '\\footskip=0.6in'
      write(lunOut,*) ' '
      write(lunOut,*) ' %%%%%% COLORS'
      write(lunOut,*) '\\newcommand{\\blue}{\\color{blue}}'
      write(lunOut,*) '\\newcommand{\\red}{\\color{red}}'
      write(lunOut,*) '\\newcommand{\\green}{\\color{green}}'
      write(lunOut,*) '\\newcommand{\\magenta}{\\color{magenta}}'
      write(lunOut,*) '\\newcommand{\\cyan}{\\color{cyan}}'
      write(lunOut,*) '\\newcommand{\\yellow}{\\color{yellow}}'
      write(lunOut,*) ''
      write(lunOut,*) '\\newcommand{\\hbrk}{\\hfill \\break}'
      write(lunOut,*) ''
      write(lunOut,*) ''
      write(lunOut,*) '\\newcommand{\\referenceA}{\\rm }'
      write(lunOut,*) '\\newcommand{\\referenceB}{\\rm  }'
      write(lunOut,*) '\\newcommand{\\referenceC}{\\rm   }'
      write(lunOut,*) ''
      write(lunOut,*) '\\pagestyle{myheadings}'
      write(lunOut,*) ''
      write(lunOut,*) ''
      write(lunOut,*) ''
      write(lunOut,*) '\\begin{document}'
      write(lunOut,*) ''
      write(lunOut,*) ''

      call somCom(
     >      nCell,aK,xi,pf,r1,r2,B2,B1,gap,p1,p2,T1,T2, 
     > Brho1, Brho2,Trev1,Trev2,ttf, Af, driftI,driftX,prtcl)
      write(6,fmt='(a,i4,1p,7g12.4)') 
     > 'nCell, aK, xi, pf, r2, gap, T1, T2 :', 
     >  nCell, aK, xi, pf, r2, gap, T1, T2

      write(lunOut,*) '  \\begin{figure}'
      write(lunOut,*) '  \\hfill F.~M\\''eot'
      write(lunOut,*) '  \\begin{center}'
      write(lunOut,*) '{\\fbox{\\bf \\Large Spiral ring}}'
      write(lunOut,*) '  \\end{center}'
      write(lunOut,*) ''
      write(lunOut,*) ' \\vspace{4ex} '
      write(lunOut,*) ''
      write(lunOut,*) '\\mbox{\\hspace{-.1\\linewidth} '
      write(lunOut,*) '\\begin{minipage}{0.6\\linewidth} '
      write(lunOut,*) 
     >   '\\includegraphics[width=12cm]{gnuplotRing.eps}'
      write(lunOut,*) '\\end{minipage}\\hspace{-2ex} '
      write(lunOut,*) '\\begin{minipage}{0.5\\linewidth} '
      write(lunOut,*) '  \\begin{center}'
c      write(lunOut,*) '    \\begin{table}'
      write(lunOut,*) '  \\large'
c      write(lunOut,*) '  \\caption{\\label{TabParam}'
      write(lunOut,*) '{\\bf Parameters, including '
      write(lunOut,*) '            sample K/$\\xi$ values}'
c      write(lunOut,*) '  }'
      write(lunOut,*) '    \\begin{tabular}{lcl}'
      write(lunOut,*) ' \\\\'
      write(lunOut,fmt='(a,1p,I2,a)') ' Nb cells &',nCell,' \\\\'
      write(lunOut,fmt='(a,1p,e11.3,a)') ' K &',aK,' \\\\'
      write(lunOut,fmt='(a,1p,e11.3,a)') ' $\\xi$ &',xi,'& (deg.)  \\\\'
      write(lunOut,fmt='(a,1p,g11.3,a)') ' pf &',pf,' \\\\'
      write(lunOut,fmt='(a,1p,e11.3,a)') ' gap &',gap,'& (m)  \\\\'
      write(lunOut,fmt='(a,3(f8.4,a))') ' $r_1/r_2/\\Delta r$&',
     >                             r1,' / ',r2,' / ',r2-r1,'& (m)  \\\\'
      write(lunOut,fmt='(a,a6,a1,1p,2(g11.3,a))') '$E_1~/~E_2$, ',prtcl,
     >                   '&',T1/1.e6,' / ',T2/1.e6,'& (MeV)  \\\\'
      write(lunOut,fmt='(a,1p,3(g12.4,a))') ' $p_1~/~p_2~/~ratio$ &',
     >               p1/1.e6,' / ',p2/1.e6,' / ',p2/p1,'& (MeV/c) \\\\'
      write(lunOut,fmt='(a,1p,2(g11.3,a))') 
     >                            ' $B\\rho_1~/~B\\rho_2$ &',
     >                       Brho1,' / ',Brho2,' &(T.m) \\\\'
      write(lunOut,fmt='(a,1p,2(g12.4,a))') 
     >                            ' $B_1~/~B_2$ &',
     >                       B1,' / ',B2,' &(T) \\\\'
      write(lunOut,fmt='(a,2(f7.2,a))') 
     >                            ' $Trev_1~/~Trev_2$ &',
     >                       Trev1*1e9,' / ',Trev2*1e9,' &(ns) \\\\'
      write(lunOut,fmt='(a,3(f9.4,a))') 
     >                            ' $Frev_1~/~Frev_2/ratio$ &',
     >    1e-6/Trev1,' / ',1e-6/Trev2,' / ',Trev1/Trev2,' &(MHz) \\\\'
      write(lunOut,fmt='(a,1p,g11.3,a)') 
     >                    ' Dip. sector angle &',Af*r2d,'& (deg.) \\\\'
      write(lunOut,fmt='(a,1p,g11.3,a)') 
     >                    ' Dip. bend angle &',ttf*r2d,' & (deg.) \\\\'
      write(lunOut,fmt='(a,1p,g11.3,a)') 
     >      ' Drift L, inj. &',driftI,'-2$\\times$0.15 & (m) \\\\'
      write(lunOut,fmt='(a,1p,g11.3,a)') 
     >      ' Drift L, xtr. &',driftX,'-2$\\times$015 & (m) \\\\'
      write(lunOut,*) '    \\end{tabular}'
c      write(lunOut,*) '    \\end{table}'
      write(lunOut,*) '  \\end{center}'
      write(lunOut,*) '\\end{minipage} '
      write(lunOut,*) '}'

      write(lunOut,*) ''
      write(lunOut,*) ' \\vspace{4ex} '
      write(lunOut,*) '  \\begin{center}'
      write(lunOut,*) '{\\bf \\large Optical functions, } \\\\'
      write(lunOut,fmt='(1p,2(a,g12.4),a)') ' K = ',aK,',   $\\xi$ = ',
     >                         xi,' deg.~:    \\\\'
      write(lunOut,*) ' \\vspace{-2ex} '
      write(lunOut,*) 
     >   '\\includegraphics[angle=-90,width=13cm]{mad.ps}'

      open(unit=lunIn,file='scanKXi.out')
      call go2end(lunIn)
      do i = 1, 2
        backspace(lunIn)  
      enddo
      read(lunIn,fmt='(a)') txt132    
      write(*,*) '  geneTex , txt132 : ', txt132    
      write(lunOut,fmt='(a)') ' '
      write(lunOut,*) ' \\vspace{-2ex} '
      write(lunOut,*) ' nb. cells,   K,   xi,   Qx,   Qy,  ',
     >  '    Maxima  BetX,  BetY,  Dx :'
      write(lunOut,fmt='(a)') ' '
       write(lunOut,fmt='(a)')txt132(debstr(txt132):finstr(txt132)-6)
      write(lunOut,*) '  \\end{center}'
      write(lunOut,fmt='(a)') ' '
      close(lunIn)
      write(lunOut,*) '  \\end{figure}'

      write(lunOut,*) '\\clearpage'

      write(lunOut,*) '  \\begin{center}'
      write(lunOut,*) '{\\bf \\large Extrema of betatron functions}\\\\'
      write(lunOut,*) '{\\large including sample working point}'
      write(lunOut,*) '  \\end{center}'
      write(lunOut,fmt='(a)') ' '
      write(lunOut,*) '\\mbox{\\hspace{-.1\\linewidth} '
      write(lunOut,*) '\\includegraphics[width=9cm]{gnuplotBtaXMax.eps}'
      write(lunOut,*) '\\includegraphics[width=9cm]{gnuplotBtaZMax.eps}'
      write(lunOut,*) '}'
      write(lunOut,*) ''
      write(lunOut,*) '\\vspace{8ex}'
      write(lunOut,*) ''
      write(lunOut,*) '  \\begin{center}'
      write(lunOut,*) '{\\bf \\large Stability domain} \\\\'
      write(lunOut,*) '{\\large including sample working point}'
      write(lunOut,*) '  \\end{center}'
      write(lunOut,fmt='(a)') ' '
      write(lunOut,fmt='(a)') '{\\bf \\large \\hfill Left : K, $\\xi$ 
     >         \\hfill      Right : $\\nu_x,~\\nu_z$ \\hfill }'
      write(lunOut,fmt='(a)') ' '
      write(lunOut,*) '\\mbox{\\hspace{-.1\\linewidth} '
      write(lunOut,*) '\\includegraphics[width=9cm]{gnuplotKxi.eps}'
      write(lunOut,*) '\\includegraphics[width=9cm]{gnuplotQxQz.eps}'
      write(lunOut,*) '}'
      write(lunOut,*) ''
      write(lunOut,*) ''
      write(lunOut,*) '\\end{document}'
 
      close(lunOut)

      call system('latex log ; dvips log')
      call system('gv --scale=1 log.ps ')
      stop
      end

      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
      implicit double precision (a-h, o-z)
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

      subroutine readat(lunIn,nCell, aK, xi, pf, r2, gap, T1, T2,prtcl)
      logical strcon
      character*132 txt132, prtcl*6
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
      read(lunIn,fmt='(a)',end=99)  txt132
      if (STRCON(txt132,'carbon',6,
     >                             IS)) then
         prtcl='carbon'
      else
         prtcl='proton'
      endif

 99   continue
      close(lunIn)
      return
      end

      subroutine somCom(
     >      nCell,aK,xi,pf,r1,r2,B2,B1,gap,p1,p2,T1,T2, 
     > Brho1, Brho2,Trev1,Trev2,ttf, Af, driftI,driftX,prtcl)
      character prtcl*(*)
      data lunIn / 7 /

      pi = 4.e0 * atan(1.e0)
      deg2rd = pi / 180.e0
      c = 2.99792458e8

      open(unit=lunIn,file='geneTex.data')

c---------------- 
      write(6,*) ' geneTex procedure now going on...' 
      call readat(lunIn,nCell, aK, xi, pf, r2, gap, T1, T2,prtcl)
      write(6,*) ' Inout data read from geneTex.data file '
c---------------- END HYPOTHESIS 
      am = 938.27231e6
      qA = 1.
      if(prtcl.eq.'carbon') then 
        am = 931.49e6
        qA = .5
      endif
      ttf = 2. * pi / nCell
      Af = pf * ttf
      p2 = sqrt(T2 * (T2 + 2.e0 *am))
      p1 = sqrt(T1 * (T1 + 2.e0 *am))
      Brho2 = p2  /c/qA
      Brho1 = p1  /c/qA
      rho2 = r2*sin(Af/2.)/sin(ttf/2.)
      r1 = r2 / (p2/p1)**(1./(1.+aK)) 
      rho1 = r1*sin(Af/2.)/sin(ttf/2.)
      B2 = Brho2 / rho2
      B1 = Brho1 / rho1
      CMag2 = 2.e0*pi * rho2
      CMag1 = 2.e0*pi * rho1
      Circ2 = CMag2 / pf
      Circ1 = CMag1 / pf
      bta2 = p2 / sqrt(p2*p2 + am*am)
      Trev2 = Circ2 / (bta2 *c)
      bta1 = p1 / sqrt(p1*p1 + am*am)
      Trev1 = Circ1 / (bta1 *c)

      flut = Circ/CMag2 - 1.e0
      rMed = (r2+r1)/2.e0
      dr = (r2-r1)/2.e0
      driftI = 2.e0*pi*(rMed-dr)*(1.e0-pf)/nCell
      driftX = 2.e0*pi*(rMed+dr)*(1.e0-pf)/nCell

      close(lunIn)
      return
      end
      subroutine go2end(lun)
      implicit double precision (a-h, o-z)
      character*132 txt
 1    continue
        read(lun,*,end=99,err=99) txt
        goto 1
 99   return
      end
