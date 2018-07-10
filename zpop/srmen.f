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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory         
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRMEN(NLOG,NL,LM,OKOPN,NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN
      CHARACTER(*) NOMFIC

C----- Synchrotron radiation

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN

      DIMENSION OX(3), WF(3), WO(6), FNR(3)
      CHARACTER(1) REP, PART
      LOGICAL EMPTY, STORE

      CHARACTER(14) GNUFIL

      SAVE GNUFIL
      SAVE KPL

      INCLUDE 'FILPLT.H'

      DATA STORE / .FALSE. /
      DATA GNUFIL / 'gnuplot.SR' /
      DATA IT /1/
      DATA KPL / 23 /

      HZ = QE / AH * 1.D3                ! Conversion KeV -> s-1

      CALL OPNDEF(NPLT,FILPLT,NL,
     >                           NOMFIC,OKOPN) 

C      IF(OKOPN) CALL XYZBR(NL,LM,KT)
      CALL READC5(KT1,KT2)
      IF(KT2.NE.KT1) THEN   
         WRITE(6,*)
     >      ' *** WARNING : SR utility will plot trajectory #KT1= ',
     >                                                KT1,' only'
         CALL READC6B(KT1,KT1)
      ENDIF
      IF(OKOPN) CALL XYZBR(NL)

      CALL SRINIT(GNUFIL,'////MENU08/16////',
     >                   PART,R0,Q,AM,OX,WF,WO,FNR,STORE)

 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,*)
      WRITE(6,*) '  Press Return   for Menu' 
      READ(5,FMT='(A1)',ERR=20) REP 

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

C      QM = Q*CL*CL/(QE*AM*1.D6)                    ! q/m, MKSA units
      FE = Q * CL * 1.D-7     ! = Q/4.Pi.Epsilon0.c;  gives  E in V/m
      YNRM = FNR(1) / FNR(2) * FNR(3)

      WRITE(6,104) NOMFIC,IT,PART
 104  FORMAT(//,3X,60('-'),//,20X,' SYNCHROTRON RADIATION MENU:',/,
     1/,5X,' 1  Open input data file - current is ',A,
     2/,5X,' 2  PARTICLE # AND TYPE ( ',
     >                          I3,', ',A1,' )',
     3/,5X,' 3  Normalization factors (I/q, etc.) ',
     5/,5X,' 5  Electric field',
     6/,5X,' 6  Spectrum',
     7/,5X,' 7  Integrals',
     8/,5X,' 8  Grahic  menu  ',
     9/,5X,' 9  Exit  this  menu ',
     2/,5X,'12  Erase  display ',
     2/,5X,'13  Analytical models ',
     >/,3X,60('-'),//)

      IF(.NOT. OKOPN) CALL OPNWRN(1)

      WRITE(6,100)
 100  FORMAT('$  YOUR  CHOICE : ')
      READ(5,108,ERR=21) IOPT
108   FORMAT(I2)

      GOTO ( 1, 2, 3,21, 5, 6, 7, 8,99,21,21,12,13) IOPT  
      GOTO 21
 
 1    CONTINUE
C----- Open data file 
      CLOSE(NL)
      NL=NPLT
      WRITE(6,*) ' Give input data file name'
      WRITE(6,*) '       - default will be :',FILPLT
      READ(5,FMT='(A)') NOMFIC
      IF( EMPTY(NOMFIC) ) NOMFIC = FILPLT
      WRITE(6,*)
      WRITE(6,*) ' Trying to open ',NOMFIC,' ...' 
      OPEN(UNIT=NL, FILE=NOMFIC,STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        WRITE(6,*) '  Cannot open ',nomfic
      ELSE
        OKOPN = .TRUE.
        WRITE(6,*) 
        WRITE(6,*) ' Wait... '
        WRITE(6,*) '   Storing trajectoriy, fields, etc., '
        WRITE(6,*) '     as being read from ',nomfic
        WRITE(6,*) 
C        CALL XYZBR(NL,LM,KT)
        CALL XYZBR(NL)
      ENDIF
C      OKECH = .FALSE.
      GOTO 20

 2    CONTINUE
 22   WRITE(6,FMT=
     >'(/,''  # of the particle to analyse ( 1-'',I6,'' ): '')') MXT
      READ(5,*,ERR=22) IT
      IF(IT .LE. 0 .OR. IT .GT. 200) GOTO 2
 23   WRITE(6,FMT='(/,''  TYPE OF PARTICLE ( E, P ): '')') 
      READ(5,FMT='(A1)',ERR=23) PART
      IF    (PART .EQ. 'P' .OR. PART .EQ. 'p') THEN
        PART = 'P'
        AM = .93827231D3
        R0 = 1.53469852D-18
        Q = 1.60217733D-19
      ELSEIF(PART .EQ. 'E' .OR. PART .EQ. 'e') THEN
        PART = 'E'
        AM =  .51099906D0
        R0 =  2.81794092D-15
        Q = 1.60217733D-19
      ELSE
        GOTO 23
      ENDIF
      GOTO 21

 3    CONTINUE
        WRITE(6,*) '  FRN(1-3) = ',(FNR(I),I=1,3)
        WRITE(6,*) '  FNR(1)/FNR(2)*FNR(3) = YNRM = ', YNRM

 31     WRITE(6,*) '  Read zpop.dat ? (n/y)'
        READ(5,FMT='(A1)',ERR=31) REP 

        IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
          CALL SRINIT(GNUFIL,'////MENU08/16////',
     >                       PART,R0,Q,AM,OX,WF,WO,FNR,STORE)
          WRITE(6,*) '  Particle type : ',PART
          WRITE(6,*) '  Value of AM, Q, R0 after read :'
          WRITE(6,*) '   R0, Q, AM = ',R0,Q,AM
          WRITE(6,*) '  Value of normalisation coeffs after read :'
          WRITE(6,*) '  FRN(1-3) = ',(FNR(I),I=1,3)
          WRITE(6,*) '  FNR(1)/FNR(2)*FNR(3) = YNRM = ',
     >      FNR(1) / FNR(2) * FNR(3)
        ENDIF
      GOTO 20

 5    CONTINUE
        CALL SREFM(NLOG,NL,LM,OX,Q,AM,FE,NOMFIC,OKOPN,*1)
      GOTO 21

 6    CONTINUE
        CALL SRDWM(NLOG,NL,LM,KPL,OX,WF,Q,AM,FE,HZ,YNRM,
     >    NOMFIC,OKECH,OKOPN,*1)
      GOTO 21
 
 7    CONTINUE
        CALL SRINM(NLOG,NL,LM,KPL,OX,WF,WO,Q,AM,FE,HZ,YNRM,
     >    NOMFIC,OKOPN,GNUFIL,STORE,*1)
      GOTO 21
 
 8    CONTINUE
        CALL MENVCF
      GOTO 21

 12   CONTINUE
        CALL CLSCR
      GOTO 21

 13   CONTINUE
        CALL SRMODL(NLOG)
      GOTO 21

 99   RETURN 
      END
