C      implicit double precision (a-h, o-z)
      data lunIn, lunOut, lunRes, lunW / 7, 8, 10, 6 /
      data npt / 100 /

C Allows change spiral orientation
      data sign /+1./  !/-1./
c dtta0 causes a global rotation of ring and COs
      data dtta0 /0./   !!! 1.570796326794896619
c draw theoretical trajectories at Cell # (1 - nCell, 99 for all)
      data kCDraw / 3 /
c draw theoretical trajectories at Cell # (1 - nCell, 99 for all)
      data kTDraw / 0 /

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
     >    nCell, AK, ZTA, pf, r0,T0,gap, akappa, A,Q, T1,T2,mCO,map,
     >    kaseV,AT,DltACN)
      nCO = mCO
      if(nCO .lt. 1 .or. nCO .gt. 1000) nCO = 7
      write(6,*) ' Inout data read from drawSpi.data file, ', 
     >'    report can be found in drawSpi.res.' 
c---------------- END HYPOTHESIS 
      ttf = 2. * pi / nCell
      Af = pf * ttf
      p0 = sqrt(T0 * (T0 + 2.* am))
      p1 = sqrt(T1 * (T1 + 2.d0 *am))
      p2 = sqrt(T2 * (T2 + 2.d0 *am))
      r2 = r0 * (p2/p0)**(1.d0/(ak+1.d0))
      Brho2 = p2  /c
      Brho1 = p1  /c
      rho2 = r2*sin(Af/2.)/sin(ttf/2.)
      r1 = r2 / (p2/p1)**(1./(1.+aK)) 
      rho1 = r1*sin(Af/2.)/sin(ttf/2.)
      B2 = Brho2 / rho2
      CMag2 = 2.d0*pi * rho2
      Circ = CMag2 / pf
      flut = Circ/CMag2 - 1.d0
      rmed = (r2+r1)/2.d0
c      rRef = rmed
      rRef = r0
      rp1 = r1 - .2  !  rayon min pole
      rp2 = r2 + .2  !  rayon max pole
      dr = (r2-r1)/2.d0
      write(6,fmt='(a,i4,1p,7g12.4)') 
     > 'nCell, aK, zta, pf, r2, gap, T1, T2 :', 
     >  nCell, aK, zta, pf, r2, gap, T1, T2

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
      write(lunW,*) ' drift_inj = ',2.d0*pi*r1*(1.d0-pf)/nCell,
     >' m,    drift_xtr = ',2.d0*pi*r2*(1.d0-pf)/nCell,' m'  
      if(lunW.eq.6) then
        lunW = lunRes
        open(unit=lunW,file='drawSpi.res')
        goto 1 
      endif
      close(lunW)

      open(unit=lunOut,file='drawSpi.out')

      b = -sign / tan(zta*pi/180.d0)
      tta1 =  log(rp1/rRef) / b 
      tta2 =  log(rp2/rRef) / b 
      drp = (rp2-rp1)/npt


C          goto 77

Compute spiral EFBs. Can be plotted using zpop, option Menu/7/20
      do j = 1, kCDraw
        i = j-1
c        plot EFBs
        tta0 = 2.d0*pi * float(i)/float(nCell)  + dtta0
        call arcSpi(rRef,b,tta0-Af/2.,rp1,rp2,npt,lunOut)
        write(lunOut,*) '%'
        call arcSpi(rRef,b,tta0+Af/2.,rp1,rp2,npt,lunOut)
        write(lunOut,*) '%'
c        circles joining EFBs
        tta = tta0 + log(rp1/rRef) / b 
        call ACrcl(0.,0.,rp1   ,tta,-Af/2.,Af/2.,npt,lunOut)
        write(lunOut,*) '%'
        tta = tta0 + log(rp2/rRef) / b 
        call ACrcl(0.,0.,rp2   ,tta,-Af/2.,Af/2.,npt,lunOut)
        write(lunOut,*) '%'
c        plot mechanical limits
        ang = .2/rRef   !  marge de 20 cm au rayon de reference (e.g., bobine, clamp)
        call arcSpi(rRef,b,tta0-Af/2.-ang,rp1,rp2,npt,lunOut)
        write(lunOut,*) '%'
        call arcSpi(rRef,b,tta0+Af/2.+ang,rp1,rp2,npt,lunOut)
        write(lunOut,*) '%'

c        write(lunOut,*) '%'
c        write(*,fmt='(i3,5g12.4,i5,a)') 
c     >  j,tta0/deg2rd,rRef,b,rp1,rp2,npt,' j,tta0,rRef,b,rp1,rp2'
      enddo
      write(lunOut,*) '%'

   
      if(kTDraw.ne.99 .and. (kTDraw.le.0 .or. kTDraw.gt.nCell)) goto 78

C--------------------------------------------------
 77   continue
      write(*,*) ' Drawing closed orbits, busy... '

      p0 = sqrt(T0 * (T0 + 2.* am))
      p1 = sqrt(T1 * (T1 + 2.* am))
      p2 = sqrt(T2 * (T2 + 2.* am))
      dr = r2 * (1.d0 - (p1/p2)**(1.d0/(ak+1.d0)))

      ddr = 0.d0
      if(nCO.gt.1) ddr = DR / float(nCO-1) 

C      do r = r1, r2,  (r2-r1)  !!/4.
      do iCO = 1, nCO
        if(nCO.gt.1) then
          r = R2-(iCO-1.)*ddr
        elseif(nCO.eq.1) then
          r =   R2-DR/2.d0
        endif
        write(lunOut,*) '%' 
c        do j = 1, kCDraw +1 ! 1 more so to close the orbits
        do j = 1, kCDraw 
          tta = 2.d0*pi * float(j-1)/float(nCell) 
     >            + log(r/rRef) / b  + dtta0
c          write(*,fmt='(i3,3g12.4,a)') 
c     >    j,cx,cy,tta/deg2rd,'   Cell # cx,cy,tta'
          rho = r * sin(Af/2.) / sin(ttf/2.)
          oop = r * cos(Af/2.) -  rho * cos(ttf/2.)
          cx = oop * cos(tta)
          cy = oop * sin(tta)
          radius = rho
          arc = ttf
          aCma = Af/2.
          aCmi = aCma - Af
C Arc de cercle Af de courbure r centre en O
          if(kTDraw.eq.99 .or. kTDraw.eq.j) then
            call ACrcl(0.,0.,r     ,tta,aCmi,aCma,npt,lunOut)
          endif
        enddo
      enddo
C      do r = r1, r2,  (r2-r1)  !!/4.
      do iCO = 1, nCO
        if(nCO.gt.1) then
          r = R2-(iCO-1.)*ddr
        elseif(nCO.eq.1) then
          r =   R2-DR/2.d0
        endif
        write(lunOut,*) '%' 
c        do j = 1, kCDraw +1 ! 1 more so to close the orbits
        do j = 1, kCDraw 
          tta = 2.d0*pi * float(j-1)/float(nCell) 
     >            + log(r/rRef) / b  + dtta0
c          write(*,fmt='(i3,3g12.4,a)') 
c     >    j,cx,cy,tta/deg2rd,'   Cell # cx,cy,tta'
          rho = r * sin(Af/2.) / sin(ttf/2.)
          oop = r * cos(Af/2.) -  rho * cos(ttf/2.)
          cx = oop * cos(tta)
          cy = oop * sin(tta)
          radius = rho
          arc = ttf
          ama = ttf/2.
          ami = ama - ttf
C Arc de cercle ttf de rayon rho centre en O' 
          if(kTDraw.eq.99 .or. kTDraw.eq.j) then
            call ACrcl(cx,cy,radius,tta,ami ,ama ,npt,lunOut)
          endif
        enddo
      enddo
C--------------------------------------------------

 78    continue
       stop

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
     >  kaseV,AT,dtta0)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include "READAT.H"
      return
      end

      subroutine ACrcl(cx,cy,radius,aRef,ami,ama,npt,lun)
      da = (ama-ami)/(npt-1.)
      do a = ami, ama, da
        x= cx + radius * cos(a+aRef)    
        y= cy + radius * sin(a+aRef)    
        write(lun,*) x, y
      enddo
      return
      end

      subroutine arcSpi(rRef,b,tta0,rp1,rp2,npt,lunOut)
      dr = (rp2 - rp1) / (npt-1.)
      do r= rp1, rp2, dr
        tta = log(r/rRef) / b + tta0
        x = r * cos(tta)
        y = r * sin(tta)
        write(lunOut,*) x, y
      enddo
      return
      end
