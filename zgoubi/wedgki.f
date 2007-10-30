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
      SUBROUTINE WEDGKI(IO,T,Z,P,WDGA,FFXT)
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
      RHO = 1/B(1,3)  * UNIT(1)
C FM Feb. 2004
C      P = P - Z*1.D-2 * TAN(WEDG * B(1,3)*1.D2)
C      P = P - Z*1.D-2 * TAN(WEDG) * (B(1,3) / UNIT(1))
C      write(*,*) P, Z*1.D-2*( -TAN(WEDG))/rho, 
C     >          (FFXT/(6.D0*RHO*COS(WEDG)))/RHO
      P = P + Z*1.D-2*( -TAN(WEDG) + FFXT / (6.D0*RHO*COS(WEDG)) ) / RHO
C It would be better to correct wedg rather than tan(wedg), it would work better 
C for small rho. 
      RETURN
      END
