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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE PMAT(A,B,C,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------
C     COMPUTE C = A * B
C     ----------------------
      DIMENSION A(II,JJ), B(JJ,KK), C(II,KK)
      DO I = 1, II
        DO K = 1, KK
          C(I,K) = 0.D0
          DO J = 1, JJ
            C(I,K) = C(I,K) + A(I,J) * B(J,K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      
