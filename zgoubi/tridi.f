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
      SUBROUTINE TRIDI(L,D,U,COE,ND,N)
C TRIDIAGONAL MATRIX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DOUBLE PRECISION L(ND),D(ND),U(ND),COE(ND)
C FM Dec 2016
      DOUBLE PRECISION L(ND),D(ND),U(ND),COE(0:ND)
      M = N - 1
      DO I = 1,M
         L(I+1) = L(I+1) / D(I)
         D(I+1) = D(I+1) - L(I+1) *U(I)
         COE(I+1) = COE(I+1) - L(I+1) *COE(I)
      ENDDO
C  THE COEFFICIENT VECTOR WILL TRANSFORM TO SOLUTION VECTOR*
      COE(N) = COE(N)/D(N)
      DO I = M,1,-1
         COE(I) = (COE(I) - U(I) *COE(I+1)) / D(I)
      ENDDO
      RETURN
      END
