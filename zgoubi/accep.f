C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE ACCEP(JJ,
     >                    EMIT,AMA,BMA,XMA,XPMA,NLIV,MXINL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C----------------------------------------------------------------------------
C     Looks for ellipse with emittance EPSPI and other parameters taken in 
C     the vicinity of the rms ellipse, that has the maximum particle content.
C----------------------------------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/UNITS/ UNIT(MXJ)
      PARAMETER (RANGE=0.2D0,ISTP=3)

      SAVE EPSPI

Compute rms ellipse
      CALL LPSFIT(JJ, 
     >                         EMIT,ALP,BET,XM,XPM)

Compute number of particles alive and number inside ellipse
      DX = XM*RANGE/ISTP
      IF(DX.EQ.0.D0) DX = 1.D-10
      DXP = XPM*RANGE/ISTP
      IF(DXP.EQ.0.D0) DXP = 1.D-10
      DA = ALP*(2.D0*RANGE)/ISTP
      IF(DA.EQ.0.D0) DA = 1.D-10
      DB = BET*(2.D0*RANGE)/ISTP
      IF(DB.EQ.0.D0) DB = 1.D-10
      MXINL = -1
      DO 10 IXX = -ISTP,ISTP
       XX = XM + IXX*DX
       DO 10 IXXP = -ISTP,ISTP
        XXP = XPM + IXXP*DXP
        DO 10 IAA = -ISTP,ISTP
         AA = ALP + IAA*DA
         DO 10 IBB = -ISTP,ISTP
          BB = BET + IBB*DB
          CALL CNTINL(JJ,EPSPI,AA,BB,XX,XXP,
     >                                         NLIV,NINL)
          IF(NINL .GT. MXINL) THEN
            MXINL = NINL   
            AMA = AA
            BMA = BB
            XMA = XX
            XPMA = XXP
          ENDIF
 10   CONTINUE
C      RATIN = DBLE(MXINL)/DBLE(IMAX)
      EMIT = EPSPI
      RETURN
      ENTRY ACCEPW(EPSPIN)
      EPSPI = EPSPIN
      RETURN
      END
