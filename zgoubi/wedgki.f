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
      SUBROUTINE WEDGKI(IO,T,Z,P,WDGA,FINT,GAP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "MAXCOO.H"
      COMMON/UNITS/ UNIT(MXJ)
      IF    (IO .EQ. 1) THEN 
C-------- Entrance
         WEDG = WDGA + T
      ELSEIF(IO .EQ. 2) THEN 
C-------- Exit
         WEDG = WDGA - T
      ENDIF
C Rho corresponds to the total field seen by the particle, i.e. contribution 
C         of all Bi components in case of combined function (e.g., 'MULTIPOL' with dip+quad+...)
      RHO = 1.D0/B(1,3) 
      PSI = FINT * GAP/RHO * (1.D0 + SIN(WEDG)**2)/COS(WEDG)
      P = P - Z * TAN(WEDG -PSI) / RHO 
CCCCCC      P = P + Z*1.D-2*( -TAN(WEDG) + FINT / (6.D0*RHO*COS(WEDG)) ) / RHO
      RETURN
      END
