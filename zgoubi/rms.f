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

      SUBROUTINE RMS(JJ,
     >                              SQ,AL,B,XM,XPM,rx)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT)
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)


C      c = 299 792 458
       C = CL

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
        rx=0.0
        NPTS = 0

        DO 21 I=1,IMAX
C---- in case a particle exits an element its IEX flag is set to -1
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
c          WRITE(*,*) F(7,I)
          X   = F(J1,I)*UNIT1
          XP  = F(J2,I)*UNIT2
          IF(J2.EQ.1) THEN
C--------- To get kineticE
              P = BORO*CL9 *XP * AMQ(2,I)
              XP= SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I)
          ENDIF

C---------  Condition to calculate the longitudinal rms beamsize
        IF (J1.EQ.7) THEN
             P = BORO*CL9 *F(1,I) * AMQ(2,I)
             ENRG=SQRT(P*P + AMQ(1,I)*AMQ(1,I))
             BTA = P/ENRG
             Fac= (BTA * c * 1e-4 * 0.5)**2  ! THE 0.5 FACTOR IS DUE TO THE DEFINITION OF rx=2*sqrt(X2)
        ELSE
             Fac = 1.0
        ENDIF
C--------------------------
          X2  = X2 + Fac * (X-YM)**2
          XP2 = XP2 + (XP-YPM)**2
          XXP = XXP + (X-YM)*(XP-YPM)
 26     CONTINUE
        X2  = X2/NPTS
        rx = 2*SQRT(X2)
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
        SQ = 4*SQRT(X2*XP2-XXP*XXP)
        IF(SQ .GT. 0.D0) THEN
          B =  X2/SQ
          AL =  -XXP/SQ
        ENDIF

      RETURN
      END
