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
      SUBROUTINE ITAVAN(FONC,N,X,XMIN,XMAX,Y,V,F0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),XMIN(*),XMAX(*),Y(*),V(*)

      EXTERNAL FONC

1     CONTINUE
      IF(IFMES(0) .EQ. 0) THEN
         DO 2 I=1,N
            Y(I)=X(I)
            IF(V(I) .NE. 0) THEN
               X(I)=X(I)+V(I)
               IF     (X(I).LE.XMIN(I)) THEN
                  V(I)=0.0
                  X(I)=XMIN(I)
               ELSE IF(X(I).GE.XMAX(I)) THEN
                  V(I)=0.0
                  X(I)=XMAX(I)
               ENDIF
            ENDIF
2        CONTINUE
         CALL CPTFCT(FONC,F1)
         IF(F1.LT.F0) THEN
            F0=F1
            GOTO 1
         ELSE
            DO 3 I=1,N
               X(I)=Y(I)
3           CONTINUE
         ENDIF
      ENDIF
      RETURN
      END
