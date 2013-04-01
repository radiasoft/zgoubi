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
      FUNCTION BESKNU(I,Q,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------------
C     I=0 : modified Bessel_k function with fractional index. 
C     I = 1 : its integral from x to infnty. 
C     After V. O. Kostroun, NIM 172 (1980) 371-374.
C     -------------------------------------------------------
      PARAMETER (AH = 0.5D0)
      PARAMETER (EPSI = 1.D-5)

      GOTO (1,2) I+1

 1    CONTINUE
      B = 0.5D0 * EXP(-X)
      R = 0.D0
      DB = 1.D10
      DO WHILE (DB .GT. EPSI)
        R = R + 1.D0
        DB = EXP( -X * COSH(R * AH)) * COSH(Q * R * AH)
        B = B + DB
      ENDDO            

      GOTO 98

 2    CONTINUE
      B = 0.5D0 * EXP(-X)
      R = 0.D0
      DB = 1.D10
      DO WHILE (DB .GT. EPSI)
        R = R + 1.D0
        DB = EXP( -X * COSH(R * AH)) * COSH(Q * R * AH)/COSH(R * AH)
        B = B + DB
      ENDDO

 98   CONTINUE
      BESKNU = B * AH 
      RETURN
      END      
