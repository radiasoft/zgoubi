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
C1       2         3         4         5         6         7
C2345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE UNITR(KXI,KYI,
     >                         UXO,UYO)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (MXJ=7)
      COMMON/UNITS/ UNIT(MXJ-1)

      IF(KXI.EQ.1) THEN
        UXO = UNIT(6)
      ELSEIF(KXI.LE.MXJ) THEN
        UXO = UNIT(KXI-1)
      ELSE
        UXO = 1.D0
      ENDIF
      IF(KYI.EQ.1) THEN
        UYO = UNIT(6)
      ELSEIF(KYI.LE.MXJ) THEN
        UYO = UNIT(KYI-1)
      ELSE
        UYO = 1.D0
      ENDIF
      RETURN
      END
