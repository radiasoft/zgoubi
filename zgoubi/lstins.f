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
C
C     F(P(1:N)) IS AN ORDERED LIST, LOW TO HIGH
C     SWAP P(I) WITH ONE OF P(1:N) SO THAT F(P(1:N)) IS ORDERED, AND
C     CONTAINS THE SMALLEST N OF {F(P(1:N)),F(P(I))}, PLACING
C     THE REMAINING P(I) INTO P(I)
C
C     F (DP(N),IN)
C     I (INT,IN)
C     P (INT(N),IN)
C     N (INT,IN)
C
      SUBROUTINE LSTINS(F,I,P,N)

      IMPLICIT NONE
      DOUBLE PRECISION F(*)
      INTEGER I,N,P(*),PT
      INTEGER J

      IF (N.EQ.0) RETURN
      IF (F(P(N)).LE.F(P(I))) RETURN
      J = N
      PT = P(N)
 1000 IF (J.LE.1) GOTO 1100
      IF (F(P(J-1)).LE.F(P(I))) GOTO 1100
      P(J) = P(J-1)
      J = J-1
      GOTO 1000

 1100 P(J) = P(I)
      P(I) = PT

      END
