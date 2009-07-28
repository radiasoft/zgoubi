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
      SUBROUTINE MINONM(N,X,P,V,Y,XI,F0,FINI)

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X(*),P(*),V(*),Y(*),XI(*),F0,FINI
      INTEGER I
      DOUBLE PRECISION NMFINI,T
      DOUBLE PRECISION NMFONC
      EXTERNAL NMFONC
      DOUBLE PRECISION NMMIN, KO
      T = NMFINI(N)
      CALL CPTINI
      FINI = 3.40282347D+38
      DO 1000 I=1,N
         V(I)=X(I)
         XI(I)=X(I)
         P(I)=1D-3*(X(I+2*N)-X(I+N))
 1000 CONTINUE
C      CALL NMMIN(NMFONC,V,P,N)
      KO = NMMIN(NMFONC,V,P,N)
      DO 2000 I=1,N
         X(I)=V(I)
 2000 CONTINUE
      F0 = NMFONC(V)

      END
