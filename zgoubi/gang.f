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
      FUNCTION GANG(IA,PX,KAXE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PX(*)

      SAVE PXY

      GANG = 0.D0
      GOTO (1,2) IA
 
 1    CONTINUE
      IF    (KAXE .EQ. 1) THEN
        PXY = SQRT(PX(1)*PX(1)+PX(2)*PX(2))
C       .... 0 < T < PI :
        GANG = ACOS(PX(1)/PXY)*1000.D0
C       .... -PI < T < PI :
        IF(PX(2) .LT. 0.D0) GANG = -GANG
      ELSEIF(KAXE .EQ. 2) THEN
C       .... 0 < T < PI :
        GANG = ACOS( PX(1)/PX(5) )*1000.D0
      ENDIF
      RETURN
 
 2    CONTINUE
      IF    (KAXE .EQ. 1) THEN
C       .... 0 < P < PI/2
        GANG = ACOS(PXY/PX(5))*1000.D0
      ELSEIF(KAXE .EQ. 2) THEN
C       .... 0 < P < PI
        GANG = 1000.D0*ACOS( PX(2)/SQRT(PX(2)*PX(2)+PX(3)*PX(3)) )
      ENDIF
      IF(PX(3) .LT. 0.D0) GANG = -GANG
      RETURN
 
      END
