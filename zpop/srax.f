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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRAX(OKECH,KSC,XMI,XMA,YMI,YMA,XDI,YDI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH
      CHARACTER XDI*(*), YDI*(*)

        WRITE(6,*) '               Xmin-max :',XMI,XMA,' ',XDI
        WRITE(6,*) '               Ymin-max :',YMI,YMA,' ',YDI

        IF(OKECH) THEN
C--------- Plot axis
          IF(KSC .EQ. 3) THEN
C----------- log-log axis. 31=without grid, 32=with grid
            MOD = 32
          ELSE
C----------- lin-lin axis. 1=without grid, 2=with grid
            MOD = 2
          ENDIF 
          CALL TRAXES(XMI,XMA,YMI,YMA,MOD)
        ENDIF
      RETURN
      END
