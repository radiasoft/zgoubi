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
      SUBROUTINE ELLIP(AK2,AKP2,C2,CP2,AK,E,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      REAL*8 K2,KP2,K
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

C      REAL*8    ALPHA,BETA,DELTA,EPSI,TILDA,SIGMA
C      REAL*8    ALPHA1,BETA1,DELTA1,EPSI1,TILDA1,SIGMA1
C      INTEGER I

      ALPHA = 1.D0
      BETA  = SQRT(AKP2)
      DELTA = CP2/SQRT(AKP2)
      EPSI  = C2/CP2
      TILDA = 0.D0
      SIGMA = 0.D0
      I     = 0

1000  IF ( ABS(DELTA - 1.D0) .GT. 1.D-4 .AND. I .LE. 25) THEN
            ALPHA1 = ALPHA
            BETA1  = BETA
            SIGMA1 = SIGMA

            ALPHA = (ALPHA1 + BETA1)/2
            BETA  = SQRT(ALPHA1 * BETA1)
            SIGMA = SIGMA + ( ((ALPHA1 - BETA1)**2)*((2.D0)**(I-1)) )

         DELTA1 = DELTA
         EPSI1  = EPSI
         TILDA1 = TILDA

         DELTA = (BETA / (4 * ALPHA))*(2 + DELTA1 + (1 / DELTA))
         EPSI  = (DELTA1 * EPSI1 + TILDA1)/(1 + DELTA1)
         TILDA = (TILDA1 + EPSI1)/2

         I = I + 1

         GOTO 1000
      ENDIF

      IF(NRES.GT.0) THEN
        IF (I .GE. 25) THEN
           WRITE (NRES,*) 'SBR ELLIP : iterations stopped at step 25'
           WRITE (6,*) 'SBR ELLIP : iterations stopped at step 25'
        ENDIF
      ENDIF

      AK = PI /(2.D0 * ALPHA)
      E = AK * (1.D0 - (AK2 + SIGMA)/2.D0)
      P = AK * (1.D0 + TILDA)

      RETURN
      END
