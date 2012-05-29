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
      SUBROUTINE EREXYZ(ER0,Z,IDE,
     >                              E,DE,DDE)
C-------------------------------------------------
C Derive derivatives d^(i+j)E/didj at (X,Y,Z) from 
C ER|z=0 and its derivatives wrt R 
C-------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER0(*)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      Z2 = Z * Z

        E(1,2) = ER0(1) + Z2*( -ER0(3)/2.D0 + Z2/24.D0*ER0(5))
        E(1,3) = Z * ( -ER0(2) + Z2/6.D0*ER0(4)) 

C        dEy/dy
        DE(2,2) = ER0(2) - 0.5D0*Z2 * ER0(4)
C        dEy/dz, dEz/dy
        DE(3,2) = Z * (-ER0(3) + Z2/6.D0*ER0(5))

C        d2Ey/dy2
        DDE(2,2,2) = ER0(3) - 0.5D0*Z2 * ER0(5)
C        d2Ey/dydz, d2Ez/dy2
        DDE(3,2,2) = -Z * ER0(4)
C        d2Ey/dz2, d2Ez/dydz
        DDE(3,3,2) = -ER0(3) + Z2/2.D0*ER0(5)

      IF(IDE .EQ. 2) RETURN

      RETURN
      END
