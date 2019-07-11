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
      SUBROUTINE CNTINL(JJ,EPSPI,ALP,BET,XM,XPM,
     >                                          LIVE,NINL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
C----- CONVERT COORD. (CM,MRD) -> (M,RD)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

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
        GAM = (1.D0+ALP*ALP) / BET

      LIVE = 0
      NINL = 0
      DO 1 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 1
          LIVE = LIVE + 1

          IF    (JJ.EQ.1) THEN
C----------- Horizontal
            Y2 = F(2,I)*UNIT(1) - XM
            T2 = F(3,I)*UNIT(2) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
          ELSEIF(JJ.EQ.2) THEN
C----------- Vertical
            Y2 = F(4,I)*UNIT(3) - XM
            T2 = F(5,I)*UNIT(4) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
          ELSEIF(JJ.EQ.3) THEN
C----------- Time-kineticE
CCCCC----------- Time-momentum
            Y2 = F(7,I) - XM
CCCCC            T2 = F(1,I)*UNIT(6) - XPM
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            T2 = SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
          ENDIF

          IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .LE. EPSPI ) NINL = NINL + 1

 1      CONTINUE
      RETURN
      END
