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
      SUBROUTINE BEAMAT(R,
     >                    F0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*),F0(6,*)
      COMMON/BEAM/ FI(6,6)

      DIMENSION FIIN(6,6)

      F0(1,1) = R(1,1)*R(1,1)*FI(1,1)-2.D0*R(1,1)*R(1,2)*FI(2,1) +
     >  R(1,2)*R(1,2)*FI(2,2)
      F0(1,2) = -(  R(1,1)*R(2,1)*FI(1,1) -
     >  ( 1.D0 + 2.D0*R(1,2)*R(2,1) )*FI(2,1) +
     >  R(1,2)*R(2,2)*FI(2,2)  )
      F0(2,1) = F0(1,2) 
      F0(2,2) = R(2,1)*R(2,1)*FI(1,1) - 2.D0*R(2,2)*R(2,1)*FI(2,1) +
     >  R(2,2)*R(2,2)*FI(2,2)
      F0(3,3) = R(3,3)*R(3,3)*FI(3,3)-2.D0*R(3,3)*R(3,4)*FI(4,3) +
     >  R(3,4)*R(3,4)*FI(4,4)
      F0(3,4) = -(  R(3,3)*R(4,3)*FI(3,3) -
     >  ( 1.D0 + 2.D0*R(3,4)*R(4,3) )*FI(4,3) +
     >  R(3,4)*R(4,4)*FI(4,4)  )
      F0(4,3) = F0(3,4) 
      F0(4,4) = R(4,3)*R(4,3)*FI(3,3) - 2.D0*R(4,4)*R(4,3)*FI(4,3) +
     >  R(4,4)*R(4,4)*FI(4,4)
      F0(1,6) = R(1,1)*FI(1,6) +R(1,2)*FI(2,6) +R(1,3)*FI(3,6) +
     >  R(1,4)*FI(4,6) +R(1,5)*FI(5,6) +R(1,6)*FI(6,6) 
      F0(2,6) = R(2,1)*FI(1,6) +R(2,2)*FI(2,6) +R(2,3)*FI(3,6) + 
     >  R(2,4)*FI(4,6) +R(2,5)*FI(5,6) +R(2,6)*FI(6,6) 
      F0(3,6) = R(3,1)*FI(1,6) +R(3,2)*FI(2,6) +R(3,3)*FI(3,6) + 
     >  R(3,4)*FI(4,6) +R(3,5)*FI(5,6) +R(3,6)*FI(6,6) 
      F0(4,6) = R(4,1)*FI(1,6) +R(4,2)*FI(2,6) +R(4,3)*FI(3,6) + 
     >  R(4,4)*FI(4,6) +R(4,5)*FI(5,6) +R(4,6)*FI(6,6) 
      RETURN

      ENTRY BEAMA1(FIIN)
C      DO 2 J=0,4,2
C        DO 2 I=1+J,2+J
C          FI(I,J+1) = FIIN(I,J+1)
C          FI(I,J+2) = FIIN(I,J+2)
C 2    CONTINUE
      DO 2 J=1,6
        DO 2 I=1,6 
          FI(I,J) = FIIN(I,J)
 2    CONTINUE
      RETURN
      END
