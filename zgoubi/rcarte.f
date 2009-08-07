C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE RCARTE(KART,IDIM,
     >                            ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     READS B-FIELD FOR CARTEMES, TOSCA, BREVOL,
C       AND E(R=0,X)-FIELD FOR ELREVOL
C     ---------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)

      CHARACTER TXT*80, STRA(2)*80

C     ... IC, IL
      READ(NDAT,*) A(NOEL,1),A(NOEL,2)

CC   rustine chicane cebaf
CC     ... BNORM...
C      READ(NDAT,*,ERR=10) A(NOEL,10), A(NOEL,11),A(NOEL,12)
C      GOTO 11
C 10   CONTINUE
C      A(NOEL,11) = 0.D0
C      A(NOEL,12) = 0.D0
C      BACKSPACE(NDAT)
C 11   CONTINUE

C----- BNORM & X-,Y-,Z-NORM
      READ(NDAT,*,ERR=8) A(NOEL,10),(A(NOEL,10+I),I=1,IDIM)
      GOTO 81
 8    CONTINUE
      WRITE(*,*) ' *** Need 4 normalisation coefficients on that line'
      CALL ENDJOB('*** Input data error in SBR RCARTE, data line #2',0)
      STOP 
 81   CONTINUE
 
C----- TITLE - Start TITLE with FLIP to get map flipped (implemented with TOSCA... to 
C                       be completed for others)
      READ(NDAT,200) TA(NOEL,1)
 200  FORMAT(A)

      IF    (IDIM .EQ. 1) THEN
        READ(NDAT,*) A(NOEL,20)
        NFIC = 1
      ELSEIF(IDIM .EQ. 2) THEN
        READ(NDAT,*) A(NOEL,20),A(NOEL,21)
        NFIC = 1
      ELSEIF(IDIM .EQ. 3) THEN
C------- TOSCA, either cartesian mesh (MOD.le.19), or polar mesh (MOD.ge.20) 
        READ(NDAT,*,ERR=50) A(NOEL,20),A(NOEL,21),KZMA,AMOD
        GOTO 51

C------- To ensure compatibility with version 3 of Zgoubi
 50     KZMA = 1
        MOD = 0
        BACKSPACE(NDAT)
 51     A(NOEL,22)=KZMA
        A(NOEL,23)=AMOD
        MOD=NINT(AMOD)
C--------------------------------------------------------
        IF    (MOD .LE. 19) THEN
C--------- Cartesian mesh
          KART = 1
          IF    (MOD .EQ. 0) THEN
C----------- 3-D map generated by  symmetrysing wrt horizontal plane
            NFIC = (KZMA/2) + 1
          ELSEIF(MOD .EQ.10) THEN
C----------- 3-D map generated by  symmetrysing wrt horizontal plane
            NFIC = (KZMA/2) + 1
          ELSEIF(MOD .EQ. 1) THEN
C----------- no symmetrization, 3-D map generated by taking 2-D maps as is
            NFIC = KZMA
          ELSEIF(MOD .EQ. 12) THEN
C--------- The all 3D map is contained in a single file
            NFIC = 1
          ENDIF

        ELSEIF(MOD .GE. 20) THEN
C--------- Cylindrical mesh. Axis is Z. The all 3D map is contained in a single file
          KART = 2
          NFIC = 1

        ENDIF
      ENDIF

C----- MAP FILE NAME(S)
      IF(NFIC.GT.1) THEN 
        DO 37 IFIC=1,NFIC
          READ(NDAT,FMT='(A)') TXT
          CALL STRGET(TXT,1,
     >                      IDUM,STRA) 
          TA(NOEL,1+IFIC) = STRA(1)
 37     CONTINUE
      ELSEIF(NFIC.EQ.1) THEN
C------- Will sum (superimpose) 1D or 2D field maps if map file name is followed by 'SUM'
        IFIC = 0
 371    CONTINUE
          STRA(2) = ' '
          IFIC = IFIC + 1
          READ(NDAT,FMT='(A)') TXT
          CALL STRGET(TXT,2,
     >                         IDUM,STRA) 
          TA(NOEL,1+IFIC) = STRA(1)
          IF(STRA(2).EQ.'SUM') THEN
            TA(NOEL,1+IFIC) = '+'//TA(NOEL,1+IFIC)(1:79)
            GOTO 371
          ENDIF
          WRITE(TA(NOEL,2+IFIC),*) 
     >                     '// Total # of maps : ',IFIC
          IF(IFIC.GT.1) 
     >          WRITE(6,*) '// Total # of maps to be summed : ',IFIC
      ENDIF

C----- DROITE(S) DE COUPURE (IA=-1, 1, 2 or 3)
      READ(NDAT,*) IA, ( A(NOEL,I),I=31,30+3*ABS(IA))

      IF(IA .GT. 3) CALL ENDJOB(
     >  '*** Error, SBR RCARTE -> input data IDRT should be < ',4)

      A(NOEL,30) = IA
C     ... IRD
      READ(NDAT,*) A(NOEL,40)
C     ... XPAS
      READ(NDAT,*) A(NOEL,50)
      ND = 50 
      IF(KART .EQ. 1) THEN
C       ... Cartesian map frame
C           KP, XCE, YCE, ALE
        READ(NDAT,*) IA,(A(NOEL,I),I=ND+10+1,ND+10+3)
        A(NOEL,ND+10) = IA
      ELSEIF(KART .EQ. 2) THEN
C       ... Polar map frame
        READ(NDAT,*) KP
        A(NOEL,ND+10) = KP
        IF( KP .EQ. 2 ) THEN
          READ(NDAT,*) (A(NOEL,I),I=ND+20,ND+20+3)
        ELSE
          READ(NDAT,*) A(NOEL,ND+20)
        ENDIF
      ENDIF

      RETURN
      END
