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
      SUBROUTINE RSPN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************
C     READS DATA FOR SPIN TRACKING
C     ****************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      INTEGER DEBSTR, FINSTR
      CHARACTER TXT20*20

C     ... INITIAL SPIN DISTRIBUTION OPTION
      READ(NDAT,*) TXT20
      READ(TXT20,*) A(NOEL,1)
      TXT20 = TXT20(DEBSTR(TXT20):FINSTR(TXT20))
      IF    (A(NOEL,1).GE.0) THEN
        READ(TXT20(1:1),FMT='(I1)') KSO
        READ(TXT20(3:20),*,END=10,ERR=10) KSO2
      ELSE
        READ(TXT20(1:2),FMT='(I2)') KSO
        READ(TXT20(4:20),*,END=10,ERR=10) KSO2
      ENDIF

 10   CONTINUE

      IF    (KSO .LE. 4) THEN
        IF     (KSO .EQ. -1) THEN
        ELSEIF(KSO .EQ. 0) THEN
        ELSEIF(KSO .EQ. 1) THEN
        ELSEIF (KSO .EQ. 4) THEN
          IF     (KSO2 .EQ. 0) THEN
            DO I=1,IMAX
              READ(NDAT,*) (SI(J,I),J=1,3)
            ENDDO
          ELSEIF(KSO2 .EQ. 1) THEN
            READ(NDAT,*) SX, SY, SZ
            DO I=1,IMAX
              SI(1,I) = SX
              SI(2,I) = SY
              SI(3,I) = SZ
            ENDDO
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
        READ(NDAT,*) A(NOEL,10), A(NOEL,11)
C       ... AL, DA = CONE ANGLE AND D-ANGLE AROUND TO, PO
        READ(NDAT,*) A(NOEL,20), A(NOEL,21)
      ENDIF 
 
      RETURN
      END
