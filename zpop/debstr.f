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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION DEBSTR(STR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER * (*) STR
C     -----------------------------------
C     Renvoie dans DEBSTR le rang du
C     premier caractere non-blanc de STR.
C     -----------------------------------
      DEBSTR=0
      LENGTH=LEN(STR)+1
1     CONTINUE
         DEBSTR=DEBSTR+1
         IF(DEBSTR.EQ. LENGTH) RETURN
         IF (STR(DEBSTR:DEBSTR).EQ. ' ') GOTO 1
      RETURN
      END
