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
      SUBROUTINE ISOFLD(OKECH,OKOPN,OKVAR,NL,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------
C     Isofield  lines,  from a map
C     read from zgoubi.map, itself
C     being field by presence of 
C     IC=2 under a field map Keyword 
C-----------------------------------
      LOGICAL OKECH, OKOPN, OKVAR
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1) 
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER REP, TXT*11, TXT132*132
      CHARACTER * 80  NOMFIC 

      DIMENSION TX(400),RY(200),HP(400,200)
      PARAMETER (MRD=9)
      DIMENSION AM(MRD),BM(MRD),CM(MRD)

      LOGICAL BINARY

      DATA NOMFIC/ 'zgoubi.map'/
      DATA REP / 'Y' /

 4    CONTINUE
      WRITE(6,*)
      WRITE(6,102) NOMFIC
 102  FORMAT('$  Read ',A,5X,
     > ' (as generated using IC=2 in zgoubi.dat) (Y/N) : ')
      READ(5,FMT='(A1)',ERR=4) REP

      IF(REP.NE. 'N' .AND. REP.NE. 'n') THEN
C        KX = 8
C        KY = 2
        KX = 48
        KY = 42
        OKVAR = .TRUE.
        IF(OKOPN) CLOSE(NL)
        OKOPN=.FALSE.
        NL=NMAP
        OPEN(UNIT=NL,FILE=NOMFIC,STATUS='OLD',ERR=97)
        OKOPN = .TRUE.
        WRITE(6,*) ' OPEN ',NOMFIC,'  OK !' 
        WRITE(6,*) ' BUSY: READING  THE  MAP  IN ',NOMFIC
        WRITE(6,*) 

C-------- Swallow 2 lines of header
        CALL HEADER(NL,2,BINARY,*98)

        READ(NL,FMT='(1X,A11,1X,I1,1P,3G15.7)',ERR=96)
     >    TXT,IDRT,(AM(I),BM(I),CM(I),I=1,IDRT)
        READ(NL,992,ERR=98) IAMAX,IRMAX ,ACN,RFR,KART,MOD
 992    FORMAT(2I3,1P,2E15.7,1X,I1,I3)
        READ(NL,FMT='(A)') TXT132 
        WRITE(6,*) ' Now reading ',TXT132 
        READ(NL,*) (TX(IA),IA=1,IAMAX)
        READ(NL,FMT='(A)') TXT132 
        WRITE(6,*) ' Now reading ',TXT132 
        READ(NL,*) (RY(IR),IR=1,IRMAX)
        READ(NL,FMT='(A)') TXT132 
        WRITE(6,*) ' Now reading ',TXT132 
        READ(NL,*) ((HP(IA,IR),IA=1,IAMAX),IR=1,IRMAX)
C        READ(NL,991) (TX(IA),IA=1,IAMAX),(RY(IR),IR=1,IRMAX)
C     >    ,((HP(IA,IR),IA=1,IAMAX),IR=1,IRMAX)
C 991    FORMAT(1P,6E18.10)

        CLOSE(NL)
        OKOPN=.FALSE.
        WRITE(6,*) 
        WRITE(6,*) ' The map is read '
        WRITE(6,*) 
        WRITE(6,*) ' Mesh is ',KPOL(KART),'  (KART =',KART,')'

C------- Convert units, from cm to m -------
        IF(KART .EQ. 1) THEN
          DO 7 IA = 1,IAMAX
            TX(IA) = TX(IA) * UNIT(1)
 7        CONTINUE
          DO 9 I=1, IDRT
            AM(I) = AM(I) / UNIT(1)
            BM(I) = BM(I) / UNIT(1)
 9        CONTINUE
        ENDIF
        RFR = RFR  * UNIT(1)
        DO 8 IR = 1,IRMAX
          RY(IR) = RY(IR) * UNIT(1)
 8      CONTINUE
C-------------------------------------------

        RFR2 = 0.D0
        CALL READCC(MOD,RFR,RFR2)

        BMAX = -1.D-10
        BMIN = 1.D10
        DO 3 IA=1,IAMAX
          DO 3 IR=1,IRMAX
            IF(HP(IA,IR) .GT. BMAX) BMAX = HP(IA,IR)
            IF(HP(IA,IR) .LT. BMIN) BMIN = HP(IA,IR)
 3      CONTINUE
        WRITE(6,*) ' IAMAX=',IAMAX,'   IRMAX=',IRMAX
     >   ,'    BMIN/MAX = ',BMIN,'/',BMAX,'  (kG)' 
        OKECH = .FALSE.
C      ELSEIF(REP.NE.'N' .AND. REP.NE.'n') THEN
C        GOTO 4
      ENDIF

      IF(.NOT.OKECH) THEN

        XMI = 1.D10
        XMA =-1.D10
        YMI = 1.D10
        YMA =-1.D10
        DO 2 IA = 1,IAMAX
          ZETA =  TX(IA)
          X = ZETA
          DO 2 IR = 1,IRMAX
            Y=RY(IR)
            IF(KART .EQ. 2) THEN
              IF    (MOD.EQ.22) THEN
C e.g., RACCAM spiral FFAG
              ELSEIF(MOD.EQ.20) THEN
C e.g., KEK 150 MeV FFAG

              ELSE
                X = Y * SIN(ZETA) 
                Y = Y * COS(ZETA) 
              ENDIF
            ENDIF
            CALL MINMAX(X,Y,XMI,XMA,YMI,YMA)
C            write(78,*) ymi,yma,ia,' isofld'
 2      CONTINUE

        WRITE(6,*) 
        WRITE(6,*) ' XMI-MAX , YMI-MAX (m) :',XMI,XMA,YMI,YMA
        IF(XMI .LT. XMA .AND. YMI .LT. YMA) THEN
          OKECH = .TRUE.
        ELSE
          WRITE(6,*) '  Error  in  scale  calculation !!'
          WRITE(6,*) '    No  magnetic  field ?'
          RETURN
        ENDIF     
C          DDX=ABS(XMA-XMI)*0.05D0
C          DDY=ABS(YMA-YMI)*0.05D0
          CALL TRAXES(XMI,XMA,YMI,YMA,2)     
C          CALL TRAXES(XMI-DDX,XMA+DDX,YMI-DDY,YMA+DDY,2)
      ENDIF

 10   CONTINUE
      WRITE(6,*)
      WRITE(6,103)
 103  FORMAT( '$ * OPTION :  NIVEAUX (1),  BREF,DBREF (2)'
     >,1X,', Change scales (4), GRAPHIC MENU (8), Quit (9): ')
      READ(5,*,ERR=10) IOPT

      IF(IOPT .EQ. 1) THEN
        WRITE(6,*)
 904    WRITE(6,104) 
 104    FORMAT('$  Number of B+dB isofield curves and range',
     >                  ' dB (+/_ , %)  :')
        READ(5,*,ERR=904) NC,DBREF
C        BREF0=BMAX/NC
        DDB = ABS(BMAX-BMIN)/NC
        BREF0=BMIN 
      ELSEIF(IOPT .EQ. 2) THEN
        NC = 1
        WRITE(6,*)
 905    WRITE(6,105)
 105    FORMAT('$  VALEURS  DE  BREF (kG)  ET  DBREF (+/- , %)  :')
        READ(5,*,ERR=905) BREF0,DBREF
      ELSEIF(IOPT .EQ. 4) THEN
 907    WRITE(6,*) ' Give Xmin-max, Ymin-max (m) :'
        READ(5,*,ERR=907) XMI,XMA,YMI,YMA
        CALL TRAXES(XMI,XMA,YMI,YMA,2) 
        GOTO 10
      ELSEIF(IOPT .EQ. 8) THEN
        CALL MENVCF
        GOTO 10
      ELSEIF(IOPT .EQ. 9) THEN
        RETURN 1
      ELSE
        GOTO 10
      ENDIF

C      CALL TXTFBG

      CALL LINTYP(1)
C      IF(IDRT.EQ. -1 .OR. IDRT.EQ. 2) THEN
      IF(IDRT.EQ. -1) THEN
        IF(AM(1).EQ. 0.D0) THEN
          X = XMI
          Y = -CM(1)/BM(1) 
          XX = XMA
          YY = Y
        ELSEIF(BM(1).EQ. 0.D0) THEN
          X = -CM(1)/AM(1)          
          Y = YMI
          XX = X
          YY = YMA          
        ELSE
          X = XMI
          Y = (-CM(1)-AM(1)*X)/BM(1)          
          XX = XMA
          YY = (-CM(1)-AM(1)*XX)/BM(1)          
        ENDIF
        CALL VECTPL(X,Y,4)        
        CALL VECTPL(XX,YY,2)        

      ELSEIF(IDRT.EQ. 1) THEN
        IF(AM(2).EQ. 0.D0) THEN
          X = XMI
          Y = -CM(2)/BM(2)          
          XX = XMA
          YY = Y
        ELSEIF(BM(2).EQ. 0.D0) THEN
          X = -CM(2)/AM(2)          
          Y = YMI
          XX = X
          YY = YMA          
        ELSE
          X = XMI
          Y = (-CM(2)-AM(2)*X)/BM(2)          
          XX = XMA
          YY = (-CM(2)-AM(2)*XX)/BM(2)          
        ENDIF
        CALL VECTPL(X,Y,4)        
        CALL VECTPL(XX,YY,2)        

      ELSEIF(IDRT.GE. 2) THEN

        DO 6 I = 1, IDRT
          IF(AM(I).EQ. 0.D0) THEN
            X = XMI
            Y = -CM(I)/BM(I)          
            XX = XMA
            YY = Y
          ELSEIF(BM(I).EQ. 0.D0) THEN
            X = -CM(I)/AM(I)          
            Y = YMI
            XX = X
            YY = YMA          
          ELSE
            X = XMI
            Y = (-CM(I)-AM(I)*X)/BM(I)          
            XX = XMA
            YY = (-CM(I)-AM(I)*XX)/BM(I)          
          ENDIF
          CALL VECTPL(X,Y,4)        
          CALL VECTPL(XX,YY,2)        
 6      CONTINUE
      ENDIF

      CALL LINTYP(9)
      BREF = BMIN
      DO 1 IC = 1,NC
C        BREF = BREF0 * IC 
        BREF = BREF + DDB
        DO 1 IA = 1,IAMAX
          ZETA =  TX(IA)
          X = ZETA
          DO 1 IR = 1,IRMAX
            Y=RY(IR)
            IF(KART .EQ. 2) THEN
              IF    (MOD.EQ.22) THEN
C e.g., RACCAM spiral FFAG
                call fbgtxt
                if(x.gt.0.2) write(89,*) x,y,' isofld'
              ELSEIF(MOD.EQ.20) THEN
C e.g., KEK 150 MeV FFAG
              ELSE
                X = Y * SIN(ZETA) 
                Y = Y * COS(ZETA) 
              ENDIF
C              X = Y * SIN(ZETA)
C              Y = Y * COS(ZETA)
            ENDIF
            IF(IA + IR .EQ. 2) CALL VECTPL(X,Y,4)
            DB = ABS((HP(IA,IR) - BREF)/BREF) * 100.D0
            IF(DB .LT. DBREF) CALL VECTPL(X,Y,2)
 1    CONTINUE
      CALL LINTYP(-1)
      CALL FBGTXT
      GOTO 10

 96   WRITE(6,*) ' *** Error during read header in ', NOMFIC
      WRITE(6,*) ' Make sure ', NOMFIC,
     >     ' has been generated using IC=2 ini zgoubi.dat'
      WRITE(6,*) ' '
      GOTO 99

 97   WRITE(6,*) ' *** Error open ',NOMFIC
      WRITE(6,*) ' '
      GOTO 99

 98   WRITE(6,*) ' *** Error during read in ', NOMFIC
      WRITE(6,*) ' Make sure ', NOMFIC,
     >     ' has been generated using IC=2 ini zgoubi.dat'
      WRITE(6,*) ' '

 99   RETURN
      END
