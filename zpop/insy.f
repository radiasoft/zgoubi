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
      SUBROUTINE INSY(
     >                IPLAN,XB0,YB0,TETA0,D,DQ,DS,DM,DSS,DPU,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ 
     >      IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      CHARACTER*10 PLAN(2)
      SAVE PLAN 

      CHARACTER*80 NOMFIC, NAMFIC
      SAVE NOMFIC
      SAVE AMAG
      LOGICAL OKOPN, CHANGE

      DATA NOMFIC / 'zgoubi.dat' /
      DATA PLAN /'Horizontal','Vertical  '/
      
      DATA IPLAN0, XB00, YB00, TETA00 / 1, 0.D0, 0.D0, 0.D0 /
      DATA AMAG / 6.D0 /
      DATA OKOPN, CHANGE / .FALSE., .FALSE./
      DATA D00,DQ00,DS00,DM00,DSS00,DPU00 /
     >   .4D0, .4D0, .2D0, .2D0, .02D0, .2D0 /

      CLOSE(NDAT)
      OPEN(UNIT=NDAT,FILE=NOMFIC,STATUS='OLD',ERR=99,IOSTAT=IOS)
      IF(IOS.NE.0) GOTO 99
      WRITE(6,FMT='(/,''  Structure data will be read from '',A)')NOMFIC
      OKOPN = .TRUE.
      REWIND(NDAT)

      IPLAN = IPLAN0
      XB0 = XB00
      YB0 = YB00
      TETA0 = TETA00
      D = D00    * AMAG
      DQ =  DQ00  * AMAG                               
      DS = DS00   * AMAG
      DM = DM00   * AMAG
      DSS = DSS00   * AMAG
      DPU = DPU00  * AMAG

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      IF(NT.EQ.-1) THEN
        CALL READC6B(1,NPTS)
      ELSE
        CALL READC6B(NT,NT)        
      ENDIF

      WRITE(6,104) NOMFIC,PLAN(IPLAN)
 104  FORMAT(5X,' MENU  -  SYNOPTIC DRAWING : ',/
     1,/,5X,' 1  Open  OPTICS  file     -  current is : ',A
     3,/,5X,' 3  Read zpop.dat for initialisations '
     4,/,5X,' 4  Plane to plot  - current is : ',A
     5,/,5X,' 5  Starting coordinates in lab. '                          
     6,/,5X,' 6  Widths and magnifying factor : '
     7,/,5X,' 7  **  Plot ** '                          
     8,/,5X,' 8  Print screen'
     9,/,5X,' 9  EXIT  THIS  MENU '
     2,/,5X,'12  ERASE  DISPLAY    '
     >,/)

      WRITE(6,100)
 100  FORMAT('$  Option  number : ')
      READ(5,108,ERR=21) IOPT
 108  FORMAT(I2)
      BINARY=BINARF(NL)
      IF(.NOT. OKOPN .AND. (IOPT.EQ.6 .OR. IOPT.EQ.7)) THEN
        CALL OPNWRN(1)           
        GOTO 20
      ELSE
        GOTO ( 1,21, 3, 4, 5, 6, 7, 8, 9,21,21,12) IOPT  
        GOTO 21
      ENDIF

 1    CONTINUE    
        CALL OPNMNL(NDAT,
     >                   NL,NOMFIC,OKOPN,CHANGE)
      GOTO 20

 3    CONTINUE
        CALL PLINIT('////MENU07////',
     >                               NL,OKOPN,CHANGE,NOMFIC)
      GOTO 20

 4    CONTINUE
      WRITE(6,FMT='('' Plane to plot is : '',A)') PLAN(IPLAN)
      WRITE(6,FMT='(''    Give new value (1/2 for H/V) : '')')
      READ(5,FMT='(I1)',ERR=4) IPLAN
      IF(IPLAN .NE. 1 .AND. IPLAN .NE. 2) GOTO 4
      IPLAN0 = IPLAN
      GOTO 21

 5    CONTINUE
      WRITE(6,FMT='(A,2G13.6)') 'Initial coordinates (m) : ',XB0,YB0
      READ(5,FMT='(2E12.4)',ERR=5) XB0,YB0
C      XB00 = XB0
C      YB00 = YB0
 51   WRITE(6,FMT='(''Initial angle (rd) : '',G13.6)') TETA0
      READ(5,FMT='(E12.4)',ERR=51) TETA0
C      TETA00 = TETA0
      WRITE(6,FMT='(''Initial angle is now  '',G13.6,'' rad'')') TETA0
      GOTO 21
 
 6    CONTINUE
      WRITE(6,161) D,DQ,DS,DM,DSS,DPU,AMAG
 161  FORMAT(5X,1P,/,'    Widths : '
     1,/,5X,' 1  Dipole     : ',G12.4
     2,/,5X,' 2  Quadrupole : ',G12.4
     3,/,5X,' 3  Sextupole  : ',G12.4
     4,/,5X,' 4  Multipole  : ',G12.4            
     5,/,5X,' 5  Drift      : ',G12.4           
     6,/,5X,' 6  Pick-up    : ',G12.4           
     >,/,5X,'             '
     7,/,5X,' 7  Magnifying factor  - current is : ', G12.4
     9,/,5X,' 9  EXIT  THIS  MENU '
     >,/)
      WRITE(6,162)
 162  FORMAT('$  Your choice : ')
      READ(5,108,ERR=6) IOPT
        GOTO ( 601,602,603,604,605,606,607, 6,21) IOPT  
        GOTO 6

 601  WRITE(6,FMT='(A,G13.6)') 'Dipole width :',D
      READ(5,*,ERR=601) D
      DOO = D
      GOTO 20
 602  WRITE(6,FMT='(A,G13.6)') 'Quadrupole width :',DQ
      READ(5,*,ERR=602) DQ
      DQOO = DQ
      GOTO 20
 603  WRITE(6,FMT='(A,G13.6)') 'Sextupole width :',DS
      READ(5,*,ERR=603) DS
      DSOO = DS
      GOTO 20
 604  WRITE(6,FMT='(A,G13.6)') 'Multipole width :',DM
      READ(5,*,ERR=604) DM
      DMOO = DM
      GOTO 20
 605  WRITE(6,FMT='(A,G13.6)') 'Box width for drift :',DSS
      READ(5,*,ERR=605) DSS
      DSSOO = DSS
      GOTO 20
 606  WRITE(6,FMT='(A,G13.6)') 'Pick-Up width :',DPU
      READ(5,*,ERR=606) DPU
      DPUOO = DPU
      GOTO 20
 607  WRITE(6,FMT='(''Geom. magnification factor : '',G13.6)') AMAG
      READ(5,*,ERR=607) AMAG
      GOTO 20

 7    CONTINUE
          RETURN
          GOTO 21

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    RETURN

 12   CONTINUE
        CALL CLSCR
      GOTO 21

 99   WRITE(6,*) 
      WRITE(6,*) ' *** Error occured upon OPEN of file ',NOMFIC
      WRITE(6,*) '   Check existence of ',NOMFIC
      RETURN 1

      ENTRY INSY1(
     >            AMAGO)
      AMAGO = AMAG
      RETURN

      ENTRY INSY2(NAMFIC)
      NOMFIC=NAMFIC
      RETURN

      ENTRY INSY4(XB0I,YB0I,TETA0I)
      XB00 = XB0I
      YB00 = YB0I
      TETA00 = TETA0I
      RETURN

      END
