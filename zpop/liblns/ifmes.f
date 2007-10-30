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
      INTEGER FUNCTION IFMES(I)
      IMPLICIT INTEGER(A-Z)
      LOGICAL ARM,ARMIT
      INTEGER EVT
      SAVE ARM,ARMIT,EVT
      DATA ARM/.FALSE./
      DATA ARMIT/.FALSE./
      DATA EVT/0/

      IFMES=0
      IF (I.LT.0.OR.I.GT.3) GOTO 10
      GOTO (1,2,3,4),I+1
C----- ARMEMENT
    1 IF (.NOT.ARM) THEN
         ARM=.TRUE.
         EVT=0
      END IF
      IF (.NOT.ARMIT) THEN
         CALL ctrlc_key(EVT)
         ARMIT=.TRUE.
      END IF
      GOTO 10
C----- DESARMEMENT
    2 CONTINUE
      ARM=.FALSE.
      EVT=0
      GOTO 10
C----- ATTENTE
    3 IF (ARM) THEN
   31    IF (EVT.NE.0) THEN
            IFMES=1
            EVT=0
         ELSE
            CALL fsleep(1)
            GOTO 31
         END IF
      END IF
      GOTO 10
C----- TEST
    4 IF (ARM) THEN
         IF (EVT.NE.0) THEN
            IFMES=1
            EVT=0
         ELSE
            IFMES=0
         END IF
      END IF
C----- SORTIE
   10 RETURN
      END    
