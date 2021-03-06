C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory       
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      OKKL=.FALSE.
      IEND=0

      IF(KL1 .GE. 0) THEN
        IF(NOEL.GE.KL1  .AND. NOEL.LE.KL2) THEN   
          OKKL=.TRUE.
        ELSE
          OKKL=.FALSE.
        ENDIF
      ELSEIF(KL1 .EQ. -1)  THEN
C--------- Take first lmnt and then every KL2 other keyword
        IF(NOEL .EQ. 1) THEN
          OKKL=.TRUE.
        ELSE
          IF((NOEL/KL2)*KL2 .EQ. NOEL) THEN
            OKKL=.TRUE.
          ELSE
            OKKL=.FALSE.
          ENDIF
        ENDIF
      ELSE
        OKKL=.FALSE.
      ENDIF
