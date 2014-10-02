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
      SUBROUTINE LIGN(K,NKAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER KBLANC,KPOINT,KOL(131)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
C      COMMON/LIGNE/ NT(200),NC(200)
      COMMON/LIGNE/ NT(MXT),NC(MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      COMMON/FAISCT/ LET(MXT)

      DATA KBLANC,KPOINT/' ','.'/

      IF(K .NE. 1)  GOTO 3
      DO 2 I=1,131
      KOL(I)=KPOINT
    2 CONTINUE
      GOTO 4
    3 CONTINUE
      DO 1 I=1,131
      KOL(I)=KBLANC
    1 CONTINUE
      DO 5 I=6,126,20
      KOL(I)=KPOINT
    5 CONTINUE
    4 CONTINUE
      IF(NKAR .EQ. 0)  GOTO 7
      DO 6 N=1,NKAR
      NK=NC(N)
      NTRA=NT(N)
      KOL(NK)=LET(NTRA)
    6 CONTINUE
    7 CONTINUE
      WRITE(NRES,100) KOL
      RETURN
  100 FORMAT(' ',129A1)
      END
