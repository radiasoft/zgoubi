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
      SUBROUTINE LPSFIT(JJ, 
     >                               SQ,A,B,XM,XPM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      COMMON/RIGID/ BORO,DPREF,DP,BR
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      COMMON/UNITS/ UNIT(MXJ)

        J1 = 2*JJ
        J2 = J1+1
        UNIT1=  UNIT(J1-1)
        UNIT2=  UNIT(J2-1)
        IF(J1.EQ.6) THEN
C--------- Time-momentum
          J1=7
          J2=1
          UNIT1=  1.D0
          UNIT2=  1.D0
        ENDIF
        XM=0.D0
        XPM=0.D0
        NPTS = 0
        DO 21 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 21
          NPTS = NPTS + 1
          X   = F(J1,I)*UNIT1
          XP  = F(J2,I)*UNIT2
          IF(J2.EQ.1) THEN
C--------- To get kineticE
              P = BORO*CL9 *XP * AMQ(2,I)
              XP= SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I)
          ENDIF
          XM  = XM + X
          XPM = XPM + XP
 21     CONTINUE
        XM = XM/NPTS
        XPM = XPM/NPTS
        YM=XM
        YPM=XPM

        J1 = 2*JJ
        J2 = J1+1
        UNIT1=  UNIT(J1-1)
        UNIT2=  UNIT(J2-1)
        IF(J1.EQ.6) THEN
C--------- Time-momentum
          J1=7
          J2=1
          UNIT1=  1.D0
          UNIT2=  1.D0
        ENDIF
        X2=0.D0
        XP2=0.D0
        XXP=0.D0
        DO 26 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 26
          X   = F(J1,I)*UNIT1
          XP  = F(J2,I)*UNIT2
          IF(J2.EQ.1) THEN
C--------- To get kineticE
              P = BORO*CL9 *XP * AMQ(2,I)
              XP= SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I)
          ENDIF
          X2  = X2 + (X-YM)**2
          XP2 = XP2 + (XP-YPM)**2
          XXP = XXP + (X-YM)*(XP-YPM)
 26     CONTINUE
        X2  = X2/NPTS
        XP2 = XP2/NPTS
        XXP = XXP/NPTS

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C gamma=d22=XP2/sqrt(DELTA),-alpha=d12=XXP/sqrt(DELTA),beta=d11=X2/sqrt(DELTA).
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
        SQ = SQRT(X2*XP2-XXP*XXP) 
        IF(SQ .GT. 0.D0) THEN
          B=  X2/SQ
          A=  -XXP/SQ
        ENDIF

      RETURN
      END
