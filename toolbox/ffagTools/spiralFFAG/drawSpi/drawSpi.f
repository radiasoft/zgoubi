C      implicit double precision (a-h, o-z)
      data lunIn, lunOut, lunRes / 7, 8, 10 /

      lunW = 6

      pi = 4.d0 * atan(1.d0)
      deg2rd = pi / 180.d0
      c = 2.99792458d8
      am = 938.27231d6

      write(6,*) 
      write(6,*) '--------------------------------------------'
      write(6,*) ' drawSpi procedure now going on...' 

c---------------- Hypothesis data : 
      call readat
     >   (lunIn,'geneMap.data', 
     >    nCell, AK, ZTA, pf, r2,T0,gap, akappa, A,Q, T1,T2,nCO,map,
     >    kaseV,AT,DltACN)

c      call readat(lunIn,nCell, aK, zta, pf, r2, gap, T1, T2)
      write(6,*) ' Inout data read from drawSpi.data file, ', 
     >'    report can be found in drawSpi.res.' 
c---------------- END HYPOTHESIS 

      write(6,fmt='(a,i4,1p,7g12.4)') 
     > 'nCell, aK, zta, pf, r2, gap, T1, T2 :', 
     >  nCell, aK, zta, pf, r2, gap, T1, T2
      ttf = 2. * pi / nCell
      Af = pf * ttf
      p2 = sqrt(T2 * (T2 + 2.d0 *am))
      p1 = sqrt(T1 * (T1 + 2.d0 *am))
      Brho2 = p2  /c
      Brho1 = p1  /c
      rho2 = r2*sin(Af/2.)/sin(ttf/2.)
      r1 = r2 / (p2/p1)**(1./(1.+aK)) 
      rho1 = r1*sin(Af/2.)/sin(ttf/2.)
      B2 = Brho2 / rho2
      CMag2 = 2.d0*pi * rho2
      Circ = CMag2 / pf
      flut = Circ/CMag2 - 1.d0
      r0 = (r2+r1)/2.d0
      rp1 = r1 - .2  !  rayon min pole
      rp2 = r2 + .2  !  rayon max pole
      dr = (r2-r1)/2.d0

 1    continue
      write(lunW,*)  
     >' drawSpi has read the following data from drawSpi.data, '
      write(lunW,*)  
     > '         nCell, K, zta, pf, r1, r2, gap, T1, T2 : '
      write(lunW,*)  ' nCell = ',nCell
      write(lunW,*)  ' K, zta = ', aK, zta ,' (deg.)'
      write(lunW,*)  ' pf = ', pf
      write(lunW,*)  ' r1, computed from r2, K, delta_p : ',r1,' m'
      write(lunW,*)  ' r2 = ',r2,' m'
      write(lunW,*)  ' gap = ',gap,' m'
      write(lunW,*)  ' T1, T2 = ',  T1, T2 ,' eV '
      write(lunW,*)  ' '
      write(lunW,*)  ' From that, '
      write(lunW,*)  ' you get rho1, rho2 = ',rho1, rho2, ' m' 
      write(lunW,*)  ' you get p1, p2 = ',p1, p2, ' MeV/c' 
      write(lunW,*)  ' you get  p2/p1= ',p2/p1 
      write(lunW,*)  ' you get Brho1, Brho2 = ',Brho1, Brho2, ' T.m' 
      write(lunW,*)  ' you get Bmax  = ',B2,' T,   at r2 =',r2,' m'
      write(lunW,*)  ' '
      write(lunW,*)  ' You get hard edge flutter = ',flut
      write(lunW,*)  ' You get theoretical nu_x / nu_z = ',
     >  sqrt(1.+aK),
     >  sqrt(-aK + flut * (1.d0+ 2.d0 * tan(zta*pi/180.d0)**2))
      write(lunW,*)  ' '
      write(lunW,*) 'It also results from that all : '
      write(lunW,*) ' drift_inj = ',2.d0*pi*(r0-dr)*(1.d0-pf)/nCell,
     >' m,    drift_xtr = ',2.d0*pi*(r0+dr)*(1.d0-pf)/nCell,' m'  
      if(lunW.eq.6) then
        lunW = lunRes
        open(unit=lunW,file='drawSpi.res')
        goto 1 
      endif
      close(lunW)

      open(unit=lunOut,file='drawSpi.out')

      b = -1.d0 / tan(zta*pi/180.d0)

      dtta0 = 0. ! pi/2. + (AT - TTF )/2.D0 - ACN 

Compute spiral EFBs. Can be plotted using zpop, option 7/20
      do j = 1, nCell
         i = j-1
        tta0 = 2.d0*pi * float(i)/float(nCell) + dtta0
C Limit angles corresponding to limit radii rp1, rp2
        tta1 =  log(rp1/r0) / b
        tta2 =  log(rp2/r0) / b
c          write(*,*) '  ttaMin, max : ', tta1, tta2
c         pause
Compute first spiral EFB
        do tta = tta1, tta2, (tta2-tta1)/100.d0
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0)  
           y = r *sin(tta + tta0)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(lunOut,*) x, y, tta, curv
        enddo
Compute mechanical limit for first  EFB
        write(lunOut,*) '  '
        dtta = -.15 / rp1     ! 10 cm at rp1 for coil, clamp...
        do tta = tta1, tta2, (tta2-tta1)/100.d0
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0 +dtta)  
           y = r *sin(tta + tta0 +dtta)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(lunOut,*) x, y, tta, curv
        enddo
Compute a circle from first to scnd spiral
        write(lunOut,*) '  '
        do tta = tta2, tta2+Af, Af/100.d0
            x = r *cos(tta + tta0)  
            y = r *sin(tta + tta0)  
            curv = 0.
            write(lunOut,*) x, y, tta, curv
        enddo
Compute scnd spiral EFB
        do tta = tta2, tta1, -(tta2-tta1)/100.d0
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0 + Af)  
           y = r *sin(tta + tta0 + Af)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(lunOut,*) x, y, tta, curv
        enddo
        ttab = tta
Compute mechanical limit for scnd  EFB
        write(lunOut,*) '  '
        dtta = .15 / rp1     ! 10 cm at rp1 for coil, clamp...
        do tta = tta2, tta1, -(tta2-tta1)/100.d0
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0 + Af + dtta)  
           y = r *sin(tta + tta0 + Af + dtta)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(lunOut,*) x, y, tta, curv
        enddo
C Get back to first point for the plotter to close the plot
        write(lunOut,*) '  '
        do tta = tta1, tta1-Af, -Af/100.d0
          x = r *cos(tta + tta0 + Af)  
          y = r *sin(tta + tta0 + Af)  
            curv = 0.
            write(lunOut,*) x, y, tta, curv
        enddo
C To avoid a line joining next pole
        write(lunOut,*) '  '
      enddo
      
C For drawing of closed orbits
      tetf = (1.d0 - pf)/pf  * Af
      
Compute closed orbits
      dtta1 = dtta0
      do r = r0-dr, r0+dr,  dr/2.
        rho = r * sin(Af/2.) / sin(ttf/2.)
        dr = r * cos(Af/2.) -  rho * cos(ttf/2.)
           write(*,*) ' dr : ', dr
        do i = 0, nCell-1
          tta1 = log(r/r0) / b + 2.d0*pi *float(i)/float(nCell) + dtta1
C          do tta = tta1, tta1+Af, Af/40.d0
          do tt = tta1, tta1+ ttf, ttf/40.d0
            x = rho *cos(tt) + dr *cos(tta1+ttf/2.)
            y = rho *sin(tt) + dr *sin(tta1+ttf/2.)
            write(lunOut,*) x, y, tta
          enddo
        enddo
Close the orbit
        tta =  log(r/r0) / b + dtta1
        x = rho *cos(tt)  + dr *cos(tta+ttf/2.)
        y = rho *sin(tt)  + dr *sin(tta+ttf/2.)
        write(lunOut,*) x, y, tta
      enddo
     
      write(lunOut,*) '  '
      write(lunOut,*) '0.   0.  ', ttaa
      write(lunOut,*) 1.2*r2*cos(tta2+Af), 1.2*r2*sin(tta2+Af)
      write(lunOut,*) '  '
      tta0 =  2.d0*pi /float(nCell) + dtta0
      write(lunOut,*) '0.   0.  ', ttaa
      write(lunOut,*) 1.2*r2*cos(ttab+tta0), 1.2*r2*sin(ttab+tta0)
      dtta =  ttab+tta0 - (tta2+Af )
      write(*,*)
      write(*,*)' Room for cavity : '
      write(*,*)' rp1*dtta = ',rp1*dtta,' m,  rp2*dtta = ',rp2*dtta,' m'

      close(unit=lunOut)

      call system('gnuplot $LOCAL/mad/structure/ffag/tools/spiralFFAG/sc
     >anKXi/gnuplotRing.cmd ')
c      call system('declare -x TOOLS="$LOCAL/mad/structure/ffag/tools" ;
c     > gnuplot $TOOLS/spiralFFAG/scanKXi/gnuplotRing.cmd ')
c      call system('gv gnuplotRing.eps &')
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

      subroutine readat(lunIn,fname, 
     >  nCell, AK, ZTA, pf, r0, T0, gap, akappa,AA,QQ,T1,T2,nCO,map, 
     >  kaseV,AT,DltACN)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include "READAT.H"
      return
      end

