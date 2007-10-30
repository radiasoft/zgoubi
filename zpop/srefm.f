C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SREFM(NL,LM,NT,IT,OX,Q,AM,FE,NOMFIC,OKOPN,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OX(*)
      CHARACTER*(*) NOMFIC
      LOGICAL OKOPN
      
C     --------------------------------------------------------
C     Plot E field from synchrotron radiation as observed in a
C     given direction OX.
C     --------------------------------------------------------
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL INPECH
      CHARACTER KVX*8,KDX*8,KVY*8,KDY*8
      CHARACTER REP

      INTEGER DEBSTR, FINSTR
      INCLUDE 'FILPLT.H'

      DATA NRMA, NRMC / 0, 0 /
      DATA IX, IY / 1, 2 /
      DATA  KVX, KDX, KVY, KDY / '  t  ', ' (s) ', '  Ey ', '(V/m)' /

C      CALL LINTYP(1)

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      WRITE(6,100) IX, IY, (OX(I),I=1,3)
 100  FORMAT(//,3X,60(1H-),//,20X,' Electric field. MENU:' ,//
     1,/,5X,' 1    X-Y coordinates to plot'
     >,/,5X,'       X: time    Obs - Part '
     >,/,5X,'                 ( 1  -   2 )                  - Now : ',I1
C     >,/,5X,'       X: time    Obs - Part -  Obs*omga_c '
C     >,/,5X,'                 ( 1  -   2  -     3 )   - Now : ',I1
     >,/,5X,'       Y: E-field    Ex,y,z '
     >,/,5X,'                    ( 1 2 3 )'
C     >,/,5X,'       Y: E-field    Ex,y,z - Ex,y,z/Emax'
C     >,/,5X,'                    ( 1 2 3 -  4 5 6 )'
     >,/,5X,'          Particle position, time  X,Y,Z, t'' '
     >,/,5X,'                                 ( 7 8 9  10) '
     >,/,5X,'                                               - Now : ',I1
     2,/,5X,' 2    Observer position w.r.t magnet frame:'
     >,/,5X,'        X, Y, Z = ',1P,3G12.4,' (M)'
     3,/,5X,' 3    Plot options'
     4,/,5X,' 4    Hand scales'
     5,/,5X,' 5    Automatic scales'
     7,/,5X,' 7    **  PLOT  Y  v.s.  X  **'
     8,/,5X,' 8    Print screen'
     9,/,5X,' 9    Exit this menu'
     X,/,5X,'10    Line style'
     2,/,5X,'12    Erase display'
     >,/,3X,60(1H-),//)

      WRITE(6,FMT='(''$  Option  number : '')')
      READ(5,FMT='(I2)',ERR=21) IOPT
      GOTO ( 1, 2, 3, 4, 5,21, 7, 8, 9,10,21,12) IOPT  
      GOTO 21

 1    CONTINUE
      IX0 = IX
      IY0 = IY
 14   WRITE(6,*) 
      WRITE(6,*) ' X variable: '
      WRITE(6,*) '    Obsrvr time , Prticl time (1-2):'
      READ(5,*,ERR=14) IX
      IF(IX .LT. 1 .OR. IX .GT. 2) GOTO 14
 13   WRITE(6,*) 
      WRITE(6,*) ' Y variable'
      WRITE(6,*) '    Electric field:     Ex, Ey, Ez '
      WRITE(6,*) '                       ( 1   2   3 )'
      WRITE(6,*) '    Particle pos., time:   X,  Y,  Z   t'''
      WRITE(6,*) '                         ( 7   8   9   10)  :'
      READ(5,*,ERR=13) IY
      IF(IY.LT. 1 .OR. IY.GT. 11) GOTO 13
C      IF(IX .NE. IX0 .OR. IY .NE. IY0 ) OKECH = .FALSE.
      IF(IX .EQ. 1) THEN
        IF(NRMA .EQ. 0) THEN
          KVX = 't'         
          KDX = '(s)'
        ELSEIF(NRMA .EQ. 1) THEN
          KVX = 'omgac*t'         
          KDX = '(rad)'
        ELSEIF(NRMA .EQ. 2) THEN
          KVX = 't'         
          KDX = '(s)'
        ENDIF
      ELSEIF(IX .EQ. 2) THEN
        KVX = 't'''         
        KDX = '(s)'
      ENDIF
      IF( IY .LE. 3 ) THEN
        IF    (IY .EQ. 1) THEN
          KVY = 'Ex'
        ELSEIF(IY .EQ. 2) THEN
          KVY = 'Ey'
        ELSE
          KVY = 'Ez'
        ENDIF
        IF(NRMA .EQ. 0) THEN
C          KVY = 'E'                 ! Ex,y or z
          KDY = '(V/m)'
        ELSEIF(NRMA .EQ. 1) THEN
          KVY = KVY(DEBSTR(KVY):FINSTR(KVY))//'/Emax' 
          KDY = '(rel.)'
        ELSEIF(NRMA .EQ. 2) THEN
          KDY = '(V/m)'
        ENDIF
      ELSE
        IF( IY .EQ. 7 ) THEN
          KVY = 'X'                
          KDY = '(m)'
        ELSEIF( IY .EQ. 8 ) THEN
          KVY = 'Y'                
          KDY = '(m)'
        ELSEIF( IY .EQ. 9 ) THEN
          KVY = 'Z'                
          KDY = '(m)'
        ELSEIF( IY .EQ. 10) THEN
          KVY = 't'''      
          KDY = '(s)'
        ELSEIF( IY .EQ. 11) THEN
          KVY = ' '                   ! Dummy ( for 1 - n.beta )
          KDY = ' '
        ENDIF
      ENDIF
      GOTO 21

 2    CONTINUE
 22   WRITE(6,*) 
      WRITE(6,*) ' Observer position  X,Y,Z (m):'
      READ(5,*,ERR=22) (OX(I),I=1,3)
C      OKECH = .FALSE.
      GOTO 21

 3    CONTINUE
        CALL PLTOPT(
     >              LM,OKECH,OKBIN,KOPT)
        IF(KOPT.EQ.8) THEN
          CALL READC5(KT1,KT2)
          IF(KT2.NE.KT1) THEN
            WRITE(6,*)
     >      ' *** WARNING : SR utility will plot trajectory #KT1= ',
     >                                                KT1,' only'
            CALL READC6B(KT1,KT1)
          ENDIF
          CALL XYZBR(NL,LM)
          GOTO 20
        ELSEIF(KOPT.EQ.88) THEN
 31       WRITE(6,*) ' Normalize Obsrvr time t to (t-t(Ey_max)/tc,'
          WRITE(6,*) '            and  Ex,y,z to Ex,y,z/Ey_max (N/Y)?'
          READ(5,FMT='(A1)',ERR=31) REP 
          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
            NRMA = 1 
 32         WRITE(6,*) ' tc and Ey_max will be computed'
            WRITE(6,*) '   once for all (1), or at each plot (2):'
            READ(5,*,ERR=32) NRMC
            IF(NRMC .NE. 1 .AND. NRMC .NE. 2) GOTO 32
          ELSE
 33         WRITE(6,*) ' Center Obsrvr time t (N/Y)'
            READ(5,FMT='(A1)',ERR=33) REP 
            IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
              NRMA = 2 
            ELSE
              NRMA=0
            ENDIF
            NRMC = 0
          ENDIF
        ENDIF
      GOTO 21
        
 4    CONTINUE
        OKECH=INPECH()
      GOTO 21
        
 5    CONTINUE
        IF(.NOT. OKOPN) CALL OPNDEF(NPLT,FILPLT,NL,
     >                                             NOMFIC,OKOPN) 
        IF(.NOT. OKOPN) RETURN 1
C        OKECH = .FALSE.
C------- Calculate scales and plot axis
        XMI = 1.D10
        XMA = -1.D10
        YMI = 1.D10
        YMA = -1.D10
        CALL SREF(1,NL,LM,NT,IT,OX,IX,IY,Q,AM,FE,NRMA,NRMC,
     >                                              GAM,R,NOC,NRD,*51)
 51     CONTINUE
        WRITE(6,*) ' CALCULATION OF SCALES FROM ',NOC,' POINTS'
        WRITE(6,*) ' Xmin-max :',XMI,XMA,' ',KDX   
        WRITE(6,*) ' Ymin-max :',YMI,YMA,' ',KDY
        IF(XMI.LT.XMA .AND. YMI.LT.YMA) THEN
C          CALL TXTFBG
          CALL CLSCR
          CALL TRAXES(XMI,XMA,YMI,YMA,-1) 
          OKECH = .TRUE.
        ELSE
          OKECH = .FALSE.
          WRITE(6,*) ' Scale min-max problem ; cannot plot!'
          IF(NOC .EQ. 0) THEN
            WRITE(6,*) 
            WRITE(6,*) '   ***  No  data  could  be  read  ***'
            WRITE(6,*) '     Check particle or element number ' 
            WRITE(6,*) '        -  Main menu, option 2  -'
          ENDIF
          RETURN 1
        ENDIF
      GOTO 20
        
 7    CONTINUE
        IF(.NOT. OKOPN) RETURN 1

        IF(.NOT. OKECH ) THEN
C--------- Calculate scales and plot axis
          CALL SREF(1,NL,LM,NT,IT,OX,IX,IY,Q,AM,FE,NRMA,NRMC,
     >                                              GAM,R,NOC,NRD,*79)
 79       CONTINUE
          WRITE(6,*)
          WRITE(6,*) ' CALCULATION OF SCALES FROM ',NOC,' POINTS'
          WRITE(6,*) ' Xmin-max :',XMI,XMA,' ',KDX
          WRITE(6,*) ' Ymin-max :',YMI,YMA,' ',KDY
          IF(XMI.LT.XMA .AND. YMI.LT.YMA) THEN
C            CALL TXTFBG
            CALL CLSCR                
            CALL TRAXES(XMI,XMA,YMI,YMA,-1) 
            OKECH = .TRUE.
          ELSE
            WRITE(6,*) ' Scale min-max problem ; cannot plot!'
            OKECH = .FALSE.
            IF(NOC .EQ. 0) THEN
              WRITE(6,*) 
              WRITE(6,*) '   ***  No  data  could  be  read  ***'
              WRITE(6,*) '     Check particle or element number ' 
              WRITE(6,*) '        -  Main menu, option 2  -'
            ENDIF
C            GOTO 20
          ENDIF
        ENDIF

       IF(OKECH) THEN
C--------- Plot radiated Electric field
          CALL SREF(2,NL,LM,NT,IT,OX,IX,IY,Q,AM,FE,NRMA,NRMC,
     >                                              GAM,R,NOC,NRD,*20)
          CALL TRKVAR(LM,NOC,KVY,KDY,KVX,KDX)
          WRITE(6,*) ' Plot  OK; END OF FILE encountered'
      ENDIF
      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    CONTINUE
C        OKECH = .FALSE.
      GOTO 99

 10   CONTINUE
        CALL TYPTRA
      GOTO 21

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 21

 99   RETURN
      END
