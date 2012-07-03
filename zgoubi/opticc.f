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
      SUBROUTINE OPTICC(LNOPTI,NOEL,KOPTIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,6),F0(6,6)

      IORD = 1
      IFOC = 0
      KWR = 0

      CALL MATRIC(IORD,IFOC,KWR)
      CALL MATRI1(
     >            R)
      CALL BEAMAT(R, 
     >              F0,PHY,PHZ)
      CALL BEAIMP(F0,PHY,PHZ)

c Store in zgoubi.OPTICS.out
      IF(KOPTIP.EQ.1) CALL OPTIMP(LNOPTI,NOEL,F0,PHY,PHZ)

      RETURN
      END
