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
      FUNCTION STRACO(NST,STRA,STR,
     >                             IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRACO
      CHARACTER STRA(*)*(*), STR*(*)     
C     -----------------------------------------------------------
C     .TRUE. if the array STRA contains the string STR
C     NST = max number of indices considered : STRA(1->NST)
C     IS = index corresponding to STR in STRA (i.e.,STRA(IS)=STR)
C     -----------------------------------------------------------
      DO 1 I = 1, NST
        IF( STRA(I) .EQ. STR ) THEN
          IS = I
          STRACO = .TRUE.
          RETURN
        ENDIF
 1    CONTINUE
      STRACO = .FALSE.
      RETURN
      END
