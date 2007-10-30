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
      SUBROUTINE MKSA(IORD,R,T,T3,T4,TX2Y,TXY2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)
      DIMENSION  T3(5,*), T4(5,*), TX2Y(5,6,*), TXY2(5,6,*)

      INCLUDE "MAXCOO.H"
      COMMON/UNITS/ UNIT(MXJ)
 
      DO 1 IA=1,6
        DO 1 IB=1,6
          R(IA,IB) = R(IA,IB)*UNIT(IA)/UNIT(IB)
C          DO 1 IC=1,IB
          DO 1 IC=1,6
            T(IA,IC,IB) = T(IA,IC,IB)*UNIT(IA)/UNIT(IC)/UNIT(IB)
 1    CONTINUE
 
      IF(IORD .EQ. 1) RETURN
 
      DO 2 IA=1,5
       DO 2 IB=1,6
        T3(IA,IB) = T3(IA,IB)*UNIT(IA)/UNIT(IB)**3
        T4(IA,IB) = T4(IA,IB)*UNIT(IA)/UNIT(IB)**4
        DO 2 IC=1,6
         TX2Y(IA,IC,IB) = TX2Y(IA,IC,IB)*UNIT(IA)/UNIT(IC)**2/UNIT(IB)
         TXY2(IA,IC,IB) = TXY2(IA,IC,IB)*UNIT(IA)/UNIT(IC)/UNIT(IB)**2
 2    CONTINUE
 
      RETURN
      END
