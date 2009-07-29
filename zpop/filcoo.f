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
      SUBROUTINE FILCOO(KPS,NOC,YZXB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YZXB(*)
      INCLUDE 'MAXNTR.H'
      PARAMETER (NTR=NTRMAX*9)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION KXYL(6)
      SAVE KXYL
      DATA KXYL / 2,3,4,5,7,20 /
C      DATA KXYL / 2,3,4,5,6,1 /
C      DATA KXYL / 2,3,4,5,18,19 /

C----------- Current coordinates
          II = 0
C----------- Initial coordinates
          IF(KPS.EQ. 0) II = 10

C if   :   DATA KXYL / 2, 3, 4, 5, 7, 20 /
C then :               x  x' z  z' t  E

          COOR(NOC,1)=YZXB(KXYL(1)+II )
          COOR(NOC,2)=YZXB(KXYL(2)+II )
          COOR(NOC,3)=YZXB(KXYL(3)+II )
          COOR(NOC,4)=YZXB(KXYL(4)+II )
          COOR(NOC,5)=YZXB(KXYL(5)+II )
          COOR(NOC,6)=YZXB(KXYL(6) )
      RETURN
      END
