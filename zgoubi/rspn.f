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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE RSPN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************
C     READS DATA FOR SPIN TRACKING
C     ****************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      INTEGER DEBSTR, FINSTR
      CHARACTER(130) TXT130
      CHARACTER(LNTA) FNAME
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE
      LOGICAL STRCON

C     ... INITIAL SPIN DISTRIBUTION OPTION
      LINE = 1
      READ(NDAT,FMT='(A)',ERR=90) TXT130
      TXT130 = TXT130(DEBSTR(TXT130):FINSTR(TXT130))
      IF(STRCON(TXT130,'!',
     >                    IS)) TXT130 = TXT130(DEBSTR(TXT130):IS-1)
      IF(STRCON(TXT130,'startAt',
     >                          IIS)) THEN
        READ(TXT130(IIS+7:FINSTR(TXT130)),*,ERR=90) A(NOEL,2)
      ELSE
        A(NOEL,2) = 0
      ENDIF
      READ(TXT130,*,ERR=90) A(NOEL,1)
      KSO2 = 0
      IF    (A(NOEL,1).GE.0) THEN
        READ(TXT130(1:1),FMT='(I1)') KSO
        IF(STRCON(TXT130,'.',
     >                      IS)) THEN
          READ(TXT130(IS+1:FINSTR(TXT130)),*,END=10,ERR=10) KSO2
        ENDIF
      ELSE
        READ(TXT130(1:2),FMT='(I2)') KSO
      ENDIF

 10   CONTINUE

      IF    (KSO .LE. 4) THEN
        IF     (KSO .EQ. -1) THEN
        ELSEIF(KSO .EQ. 0) THEN
        ELSEIF(KSO .EQ. 1) THEN
        ELSEIF (KSO .EQ. 4) THEN
          IF     (KSO2 .EQ. 0) THEN
            DO I=1,IMAX
              LINE = LINE + 1
              READ(NDAT,*,ERR=90) (SI(J,I),J=1,3)
            ENDDO
          ELSEIF(KSO2 .EQ. 1) THEN
            LINE = LINE + 1
            READ(NDAT,*,ERR=90) SX, SY, SZ
            DO I=1,IMAX
              SI(1,I) = SX
              SI(2,I) = SY
              SI(3,I) = SZ
            ENDDO
          ELSEIF(KSO2 .EQ. 2) THEN
            LINE = LINE + 1
            READ(NDAT,*,ERR=90) FNAME
            TA(NOEL,1) = FNAME
          ENDIF
        ENDIF

        IA = 0
        DO I=1,IMAX
          IA = IA + 10
          IF(IA.LE.MXD/10) THEN
            A(NOEL,IA)   = SI(1,I)
            A(NOEL,IA+1) = SI(2,I)
            A(NOEL,IA+2) = SI(3,I)
          ELSE
            IA = IA-10
            GOTO 11
          ENDIF
        ENDDO
 11     CONTINUE

        A(NOEL,9) = IA

      ELSEIF(KSO .EQ. 5) THEN
C        ... TO, PO = MEAN INITIAL PRCESSION DIRECTION
        LINE = LINE + 1
        READ(NDAT,*,ERR=90) A(NOEL,10), A(NOEL,11)
C       ... AL, DA = CONE ANGLE AND D-ANGLE AROUND TO, PO
        LINE = LINE + 1
        READ(NDAT,*,ERR=90) A(NOEL,20), A(NOEL,21)
C       ... IR random seed
        LINE = LINE + 1
        READ(NDAT,*,ERR=90) A(NOEL,30)
      ENDIF 
 
      RETURN

 90   CONTINUE
      CALL ZGKLEY( 
     >            KLE)
      CALL ENDJOB('*** Pgm rspn, keyword '//KLE//' : '// 
     >'input data error, at line #',LINE)
      RETURN
      END
