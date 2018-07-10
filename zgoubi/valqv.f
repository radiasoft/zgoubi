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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION VALQV(TIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FIRST, IDLUNI
      SAVE FIRST
      CHARACTER(10) OLDFIL
      PARAMETER (MXDAT=30000)
      DIMENSION QV(MXDAT), TIME(MXDAT)
      SAVE QV, TIME
      SAVE LNR, ICUR, TIM1

      DATA FIRST / .TRUE. /
      DATA OLDFIL / 'valqv.data' /

      IF(FIRST) THEN
        FIRST = .FALSE.
        TIM1=0.D0
        QV1=0.D0
        VALQV =  0.D0 
        ICUR=0
        IF(IDLUNI(
     >            LNR)) THEN
          OPEN( UNIT=LNR, FILE=OLDFIL, ERR=96)
        ELSE
         GOTO 95
        ENDIF
        I=0
 1      CONTINUE
          I = I+1
          IF(I.GT.MXDAT) GOTO 97
          READ(LNR,*,ERR=9,END=10) TIME(I), QV(I)
C--------- QV is transformed to Volts
          GOTO 1
      ENDIF

 9    CONTINUE

      IF(I.EQ.1) GOTO 98

 10   CONTINUE
 2    CONTINUE
        ICUR=ICUR+1
        TIM2=TIME(ICUR)
        QV2 = QV(ICUR)
        IF(TIN.GT.TIM2) THEN
          TIM1=TIM2
          QV1 = QV2
          GOTO 2
        ENDIF
      VALQV=QV1 + (QV2-QV1)/(TIM2-TIM1)*(TIN-TIM1)
      ICUR=ICUR-1
C------ Conversion from MVolts to Volts
      VALQV=VALQV*1.D6
      RETURN
 95   CALL ENDJOB('ERROR : no free unit # for '//OLDFIL,-99)
      RETURN
 96   CALL ENDJOB('ERROR upon open  old  file '//OLDFIL,-99)
      RETURN
 97   CALL ENDJOB('ERROR : max number of QV data is ',MXDAT)
      RETURN
 98   CALL ENDJOB('ERROR upon read from file '//OLDFIL,-99)
      VALQV = 0D0
      RETURN
      END
