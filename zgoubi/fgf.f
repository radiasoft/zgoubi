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
      FUNCTION  FGF(I,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(I .EQ. 2)GOTO 2
      IF(I .EQ. 3)GOTO 3
      IF(I .EQ. 4)GOTO 4
      IF(I .EQ. 5)GOTO 5
      IF(I .EQ. 6)GOTO 6
      IF(I .EQ. 7)GOTO 7
      IF(I .EQ. 8)GOTO 8
      IF(I .EQ. 9)GOTO 9
      IF(I .EQ. 10)GOTO 10
      IF(I .EQ. 11)GOTO 11
      IF(I .EQ. 12)GOTO 12
      IF(I .EQ. 13)GOTO 13
      IF(I .EQ. 14)GOTO 14
      IF(I .EQ. 15)GOTO 15
      FGF=1D0
      RETURN
2     CONTINUE
      FGF=X
      RETURN
3     CONTINUE
      FGF=Y
      RETURN
4     CONTINUE
      FGF=X**2
      RETURN
5     CONTINUE
      FGF=X*Y
      RETURN
6     CONTINUE
      FGF=Y*Y
      RETURN
7     CONTINUE
      FGF=X**3
      RETURN
8     CONTINUE
      FGF=(X**2)*Y
      RETURN
9     CONTINUE
      FGF=Y*Y*X
      RETURN
10    CONTINUE
      FGF=Y**3
      RETURN
11    CONTINUE
      FGF=X**4
      RETURN
12    CONTINUE
      FGF=(X**3)*Y
      RETURN
13    CONTINUE
      FGF=(X**2)*Y**2
      RETURN
14    CONTINUE
      FGF=X*Y**3
      RETURN
15    CONTINUE
      FGF=Y**4
      RETURN
      END
