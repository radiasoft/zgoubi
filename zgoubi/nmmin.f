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
      DOUBLE PRECISION FUNCTION NMMIN(F,X,DX,N,IARR)
C     MINIMIZE F(X), LEAVING X THAT MINIMIZES IN X
C     F (DP,IN): The function to be minimized, takes one argument X(N)
C     X (DP(N),IN/OUT): Arguments to F, initial conditions on input,
C                       arguments for minimum on output
C     DX (DP(N),IN): Initial deviations in each coordinates
C     N (INT,IN): The dimension of X
      IMPLICIT NONE
      DOUBLE PRECISION F
      EXTERNAL F
      INTEGER N
      DOUBLE PRECISION X(N),DX(N)

      INTEGER NMAX
      PARAMETER (NMAX=100)
      DOUBLE PRECISION NMMIN1,S(NMAX),Y0,Y1
      INTEGER I
      INTEGER IARR

      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.220703125D-4)

      Y1 = 3.40282347D+38
      iarr = 0

      DO 1000 I=1,N
         S(I) = X(I)
 1000 CONTINUE

 2000 Y0 = Y1
      DO 2100 I=1,N
         X(I) = S(I)
         IF (ABS(DX(I)).LT.EPS*ABS(X(I))) DX(I)=EPS*ABS(X(I))
 2100 CONTINUE
      Y1 = NMMIN1(F,S,DX,N,IARR)

      IF(IARR.EQ.1) RETURN

      IF (Y1.LT.Y0) GOTO 2000
      NMMIN=Y0

      RETURN
      END
