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
      SUBROUTINE SRDWPL(NLOG,KSC,KPL,J1,J2,J3,MX,DW,YNRM,
     >                                                         SYDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MSAM=1000)
      DIMENSION DW(MSAM,*)
      DIMENSION SYDX(*)

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB        

C      DOUBLE PRECISION SYDX, YDX

      CALL LINTYP(-1)
      HB = AH / DPI

C------- Plot curve
        CALL FBGTXT
        WRITE(6,*) ' Busy... plotting'
C        CALL TXTFBG

        DO 2 IS = 2,4
 2        SYDX(IS) = 0.D0

C------- KPL = 2,3,23 for plot of sigma, pi, or both
        IF(KPL .EQ. 23) THEN
C--------- Plot sigma & pi
          IS1 = J2
CCC To plot sigma, pi & sigma+pi :          IS2 = J3+1
          IS2 = J3
        ELSE
C--------- Plot sigma or pi
          IS1 = J1 + (KPL+1)/2
          IS2 = IS1
C          IS2 = J1 + (KPL+1)/2
        ENDIF

        DO 3 IS = IS1, IS2
          IX = 1
          X = DW(IX,J1)                          ! omga (keV)
          IF(IS .LE. J3) THEN
            Y = DW(IX,IS)                       ! dW/dNu (J/Hz)
          ELSE
            Y = DW(IX,J2) + DW(IX,J3)
          ENDIF

          IF(KSC .EQ. 2) THEN
C----------- Plot dW/dWaveL (J/m) v.s. WaveL(micro-m)
            X = 1.D-3 * AH / QE * CL / X         ! WaveL, m
            Y = Y * CL / X / X
            X = X * 1.D6                         ! WaveL, Micro-m
          ELSEIF(KSC .EQ. 3) THEN
            Y = Y / HB
          ENDIF        

C--------- Normalization, normally by I(A)/q(C) 
C          to get the power (W)  emitted by a circulating beam, 
C          from the energy radiated by 1 particle, 1 pass.
          Y = Y * YNRM

CCCCCCCCCCCCCCCCCC avril 99
C          YDX = 0.D0
          YDX = Y * (DW(2,J1) - X)*0.5D0

            ISS = IS
            IF(IS .EQ. 5 .OR. IS.EQ.8) THEN
              ISS = 2
            ELSEIF(IS .EQ. 6 .OR. IS.EQ.9) THEN
              ISS = 3
            ELSEIF(IS .EQ. 7 .OR. IS.EQ.10) THEN
              ISS = 4
            ENDIF
          SYDX(ISS) = SYDX(ISS) + YDX
          X0 = X

          CALL VECTPL(X,Y,4)
          IF(LIS .EQ. 2) CALL IMPV(NLOG,IX,X,Y,YDX,SYDX(IS),IDUM) 

          DO 31 IX = 2, MX-1
            X = DW(IX,J1)
            IF(IS .LE. J3) THEN
              Y = DW(IX,IS)
            ELSE
              Y = DW(IX,J2) + DW(IX,J3)
            ENDIF
            IF(KSC .EQ. 2) THEN
C------------- Plot dW/dWaveL (J/m) v.s. WaveL(micro-m)
              X = 1.D-3 * AH / QE * CL / X 
              Y = Y * CL / X / X
              X = X * 1.D6
            ELSEIF(KSC .EQ. 3) THEN
              Y = Y / HB 
            ENDIF        

            Y = Y * YNRM

            YDX = (X-X0) * Y

            ISS = IS
            IF(IS .EQ. 5 .OR. IS.EQ.8) THEN
              ISS = 2
            ELSEIF(IS .EQ. 6 .OR. IS.EQ.9) THEN
              ISS = 3
            ELSEIF(IS .EQ. 7 .OR. IS.EQ.10) THEN
              ISS = 4
            ENDIF
            SYDX(ISS) = SYDX(ISS) + YDX
            X0 = X

            CALL VECTPL(X,Y,2)
            IF(LIS .EQ. 2) 
     >         CALL IMPV(NLOG,IX,X,Y,YDX,SYDX(ISS),IDUM) 

 31       CONTINUE

            IX = MX
            X = DW(IX,J1)
            IF(IS .LE. J3) THEN
              Y = DW(IX,IS)
            ELSE
              Y = DW(IX,J2) + DW(IX,J3)
            ENDIF
            IF(KSC .EQ. 2) THEN
C------------- Plot dW/dWaveL (J/m) v.s. WaveL(micro-m)
              X = 1.D-3 * AH / QE * CL / X 
              Y = Y * CL / X / X
              X = X * 1.D6
            ELSEIF(KSC .EQ. 3) THEN
              Y = Y / HB 
            ENDIF        

            Y = Y * YNRM

            YDX = ((X-X0) * Y) * 0.5D0

            ISS = IS
            IF(IS .EQ. 5 .OR. IS.EQ.8) THEN
              ISS = 2
            ELSEIF(IS .EQ. 6 .OR. IS.EQ.9) THEN
              ISS = 3
            ELSEIF(IS .EQ. 7 .OR. IS.EQ.10) THEN
              ISS = 4
            ENDIF
            SYDX(ISS) = SYDX(ISS) + YDX

            CALL VECTPL(X,Y,2)
            IF(LIS .EQ. 2)  
     >        CALL IMPV(NLOG,IX,X,Y,YDX,SYDX(ISS),IDUM) 
             
 3      CONTINUE

C----- Store integral of (sigma+pi) though sigma+pi is not plotted, 
C          for printing at bottom of graphic by SRPLI
      IF(IS2.EQ.IS1+1) SYDX(4) = SYDX(2)+SYDX(3)

      RETURN
      END
