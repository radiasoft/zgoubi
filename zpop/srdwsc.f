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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRDWSC(KSC,KPL,J1,J2,J3,MX,DW,YNRM,
     >                                              OKECH,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MSAM=1000)
      DIMENSION DW(MSAM,*)
      LOGICAL OKECH

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      HB = AH / DPI

      IER = 0
      YMI = 1.D12
      YMA = -1.D12
      DO 1 IX=1,MX
        Y2 = DW(IX,J2)
        Y3 = DW(IX,J3)
        IF( KPL .EQ. 2) THEN
C--------- Plot dWy
          IF(Y2 .LT. YMI) THEN
            YMI = Y2
            XYMI = DW(IX,J1)
          ENDIF
          IF(Y2 .GT. YMA) THEN
            YMA = Y2
            XYMA = DW(IX,J1)
          ENDIF
        ELSEIF( KPL .EQ. 3) THEN
C--------- Plot dWz
          IF(Y3 .LT. YMI) THEN
            YMI = Y3
            XYMI = DW(IX,J1)
          ENDIF
          IF(Y3 .GT. YMA) THEN
            YMA = Y3
            XYMA = DW(IX,J1)
          ENDIF
        ELSEIF(KPL .EQ. 23) THEN
C--------- Plot dWy, dWz and their sum
          IF(Y2 .LT. YMI) THEN
            YMI = Y2
            XYMI = DW(IX,J1)
          ENDIF
          IF(Y3 .LT. YMI) THEN
            YMI = Y3
            XYMI = DW(IX,J1)
          ENDIF
          IF(Y2 + Y3 .GT. YMA) THEN
            YMA = Y2 + Y3
            XYMA = DW(IX,J1)
          ENDIF
        ENDIF
 1    CONTINUE
      XMI = DW(1,J1)
      XMA = DW(MX,J1)

      IF(KSC .EQ. 2) THEN
C------- Plot dW/dWaveL (MKSA) v.s. Wavelength(micro-m)
C        Conversion from omga(keV) to WaveL(micro-m)
C           Linear scales
        TEMP = 1.D-3 * AH / QE * CL / XMI * 1.D6
        XMI = 1.D-3 * AH / QE * CL / XMA * 1.D6
        XMA = TEMP
C        Conversion from dW/dNu (MKSA) to dW/dWaveL (MKSA)
        IF(XYMI*XYMA .EQ. 0.D0) THEN
          IER = -1
          OKECH = .FALSE.
          RETURN
        ENDIF
        WL =  1.D-3 * AH / QE * CL / XYMA            ! WaveL at YMA, meters
        YMA = YMA * CL / WL / WL                     ! Max(dW/dWaveL)
        WL =  1.D-3 * AH / QE * CL / XYMI
        YMI = YMI * CL / WL / WL
        IF(YMI .GT. YMA) THEN
          TEMP = YMI
          YMI = YMA
          YMA = TEMP
        ENDIF
      ELSEIF(KSC .EQ. 3) THEN
C------ Plot dN_dot / (domga/omga).dOmga  v.s. omga
C          Log scales
        YMI = YMI   / HB
        YMA = YMA   / HB
        IF(YMA .LE. 0.D0) THEN
          IER = -1
          OKECH = .FALSE.
          RETURN
        ELSE
          IF (YMI .LE. 0.D0) YMI = YMA*1.D-10 
        ENDIF
      ENDIF

C--------- Normalization, normally by I(A)/q(C) 
C          to get the power (W)  emitted by a circulating beam, 
C          from the energy radiated by 1 particle, 1 pass.
      YMI = YMI * YNRM
      YMA = YMA * YNRM
      OKECH = .TRUE.

      RETURN
      END
