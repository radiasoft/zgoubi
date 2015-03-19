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
C  -------
      SUBROUTINE REMMA(KART,IDIM,
     >                            ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------
C     READS B-FIELD FOR EMMA FFAG,
C     ---------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(80) TXT, STRA(2)

C----- IC, IL
      READ(NDAT,*) A(NOEL,1),A(NOEL,2)

C----- BNORM & X-,Y-,Z-NORM
      READ(NDAT,*,ERR=8) A(NOEL,10),(A(NOEL,10+I),I=1,IDIM)
      GOTO 81
 8    CONTINUE
      WRITE(6,*) ' *** Need 4 normalisation coefficients on that line'
      CALL ENDJOB('*** Input data error in SBR RCARTE, data line #2',0)
      STOP 
 81   CONTINUE
 
C----- TITLE - Start TITLE with FLIP to get map flipped 
      READ(NDAT,200) TA(NOEL,1)
 200  FORMAT(A)

C----- Either cartesian mesh (MOD.le.19), or polar mesh (MOD.ge.20) 
      READ(NDAT,*,ERR=50) A(NOEL,20),A(NOEL,21),KZMA,A(NOEL,23)
C      write(*,*) 'remma ',A(NOEL,20),A(NOEL,21),KZMA,A(NOEL,23)
      GOTO 51

C----- To ensure compatibility with version 3 of Zgoubi
 50   KZMA = 1
      MOD = 0
      BACKSPACE(NDAT)
 51   A(NOEL,22)=KZMA
      MOD=NINT(A(NOEL,23))
C------------------------------------------------------
      NFIC = 2
      IF    (MOD .LE. 19) THEN
C------- Cartesian mesh
        KART = 1
      ELSEIF(MOD .GE. 20) THEN
C------- Cylindrical mesh. Axis is Z. 
C MOD=22 :  a pair of data files, one for Foff one for Doff
C MOD=24 :  NN pairs of field maps, each pair corresponding to a different value 
C of the distance (d_1 ... d_NN) between quads' axis. 
        KART = 2
C        MOD=24 : Read name of the file that contains names of all map files
        IF(MOD.EQ.24) NFIC = 1
      ENDIF

      IF(MOD.EQ.0) THEN
C------- AF, AD,  DISTance between axis of quads
        READ(NDAT,*) A(NOEL,30),A(NOEL,31),A(NOEL,32)
      ELSEIF(MOD.EQ.1) THEN
C------- AF, AD, DIST1, DIST2
        READ(NDAT,*) A(NOEL,30),A(NOEL,31),A(NOEL,32),A(NOEL,33)
      ELSEIF(MOD.EQ.22) THEN
C------- AF, AD, DIST
        READ(NDAT,*) A(NOEL,30),A(NOEL,31),A(NOEL,32)
      ELSEIF(MOD.EQ.24) THEN
C------- AF, AD, DIST
        READ(NDAT,*) A(NOEL,30),A(NOEL,31),A(NOEL,32)
      ELSE
C        CALL ENDJOB('*** Input data error -> No such option MOD = ',MOD)
      ENDIF

C----- MAP FILE NAME(S)
      DO 37 IFIC=1,NFIC
        READ(NDAT,FMT='(A)') TXT
        CALL STRGET(TXT,1,
     >                       IDUM,STRA) 
        TA(NOEL,1+IFIC) = STRA(1)
C               write(*,*) ' sbr remma : ', TA(NOEL,1+IFIC) 
C               write(*,*) ' sbr remma : ', noel, ific
 37   CONTINUE

C----- DROITE(S) DE COUPURE (IA=0,-1, 1, 2 or 3)
C      READ(NDAT,*) IA, ( A(NOEL,I),I=41,43)
      READ(NDAT,*) IA, ( A(NOEL,I),I=41,40+3*ABS(IA))
      IF(IA .GT. 3) CALL ENDJOB(
     >  '*** Error, SBR RCARTE -> input data IDRT should be < ',4)
      A(NOEL,40) = IA

C     ... IRD
      READ(NDAT,*) A(NOEL,50)
C     ... XPAS
      READ(NDAT,*) A(NOEL,60)
      ND = 60 
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
