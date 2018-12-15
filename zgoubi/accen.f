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
      SUBROUTINE ACCEN(JJ,
     >                    RATIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C--------------------------------------------------------------------------
C     Computes the rms ellipse with emittance AY matched to a particle set,
C     and gives particle rate within that ellipse
C--------------------------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      SAVE EPSPI

Compute rms ellipse
      CALL LPSFIT(JJ,
     >               EMIT,ALP,BET,XM,XPM)
Compute number of particles alive and number inside ellipse
      CALL CNTINL(JJ,EPSPI,ALP,BET,XM,XPM,
     >                                    NLIV,NINL)
C        GAM = (1.D0+ALP*ALP) / BET
CC------- JJ = 1, 2 or 3
C        J1=2*JJ
C        J2=J1+1
C        IF(J1.EQ.6) THEN
C          J1=1
C          J2=6
C        ENDIF
C        NINW=0
C        NOFF = 0
C        DO 1 I=1,IMAX
C            IF(IEX(I) .LT. -1) THEN
C              NOFF = NOFF+1
C              GOTO 1
C            ENDIF
C            Y2 = F(J1,I)*UNIT(J1-1)
C            T2 = F(J2,I)*UNIT(J2-1)
C            YT = Y2*T2
C            Y2 = Y2*Y2
C            T2 = T2*T2
C            IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .LE. EPSPI) NINW=NINW+1
C 1      CONTINUE
C----- RATIN is of course at maximum 1.
C      RATIN = DBLE(NINW)/DBLE(IMAX-NOFF)
C      RATIN = DBLE(NINL)/DBLE(NLIV)
      RATIN = DBLE(NINL)/DBLE(IMAX)
      RETURN
      ENTRY ACCENW(EPSPIN)
      EPSPI = EPSPIN
      RETURN
      END
