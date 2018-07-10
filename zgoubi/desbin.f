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
      SUBROUTINE DESBIN(AM ,AM1,AM2,PX1,PX2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...........DESINTEGRATION EN 2 PARTICULES DE MASSE QCONQUE
      DIMENSION PX1(*),PX2(*)
      DATA DPI /6.28318530718D0/
 
      PX1(4)=(AM+(AM1*AM1-AM2*AM2)/AM)/2.D0
      PX2(4)=AM-PX1(4)
      PX1(5)=SQRT(PX1(4)*PX1(4)-AM1*AM1)
      PX2(5)=PX1(5)
      U=RNDM()
      COSTET=1.D0-2.D0*U
      SINTET=SQRT(1.D0-COSTET*COSTET)
      U=RNDM()
      PHI=DPI*U
      PX1(1)=PX1(5)*COSTET
      AUX=PX1(5)*SINTET
      PX1(2)=AUX*COS(PHI)
      PX1(3)=AUX*SIN(PHI)
      PX2(1)=-PX1(1)
      PX2(2)=-PX1(2)
      PX2(3)=-PX1(3)
        RETURN
          END
