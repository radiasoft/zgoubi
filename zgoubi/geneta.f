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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE GENETA(AK,E3CM,E4CM,P3,P4,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P3(*),P4(*),B(*)
C      DATA PI,DPI,RAD/3.14159265359D0,6.28318530718D0,
C     >1.745329252D-2/
C      DATA PI,DPI,RAD/3.14159265359,6.28318530718,1.745329252E-2/
 
      DPI = 8.D0 * ATAN(1.D0)
C .......3He et eta (pi) dans le CM
      U=RNDM()
      COSTET=1.D0-2.D0*U
      SINTET=SQRT(1.D0-COSTET*COSTET)
      U=RNDM()
      PHI=DPI*U
      P3(1)=AK*COSTET
      AUX=AK*SINTET
      P3(2)=AUX*COS(PHI)
      P3(3)=AUX*SIN(PHI)
      P3(4)=E3CM
      P3(5)=AK
      P4(4)=E4CM
      P4(5)=AK
      DO 20 I=1,3
        P4(I)=-P3(I)
 20     B(I)=P4(I)/P4(4)
      RETURN
       END
