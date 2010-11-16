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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION GAMMLN(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COF(6)
C      REAL*8 COF,STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF, STP / 76.18009173D0, -86.50532033D0, 24.01409822D0, 
     > -1.231739516D0, .120858003D-2, -.536382D-5, 2.50662827465D0 /
      DATA HALF, ONE, FPF / .5D0, 1.D0, 5.5D0 /
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
 11   CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END
