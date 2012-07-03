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
      SUBROUTINE OPTIMP(LUN,NOEL,F0,PHY,PHZ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F0(6,6)
      PARAMETER (PI = 4.D0*ATAN(1.D0))

      CALL SCUMR(
     >            XL,SCUM,TCUM)

 104  FORMAT(1P,13(E13.5,1X),1X,I5)
      WRITE(LUN,104) F0(1,2), F0(1,1), F0(3,4), F0(3,3) 
     >, F0(5,6), F0(5,5)
     >, F0(1,6), F0(2,6), F0(3,6), F0(4,6)
     >, PHY/(2.D0*PI), PHZ/(2.D0*PI),SCUM,NOEL
      RETURN
      END
