C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Meot
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
C  François Meot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE GETDST(NOMFIC,NBDST,
     >                               DISTL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) NOMFIC(*)
      DIMENSION DISTL(*)
      INTEGER DEBSTR, FINSTR
      CHARACTER(3) TXT3 

      DO J = 1, NBDST

        JJ = 2*J-1
        IDB = DEBSTR(NOMFIC(JJ))
        IFI = FINSTR(NOMFIC(JJ))

        I = IDB
 1      CONTINUE
          IF(NOMFIC(JJ)(I:I+2).EQ.'Dax') THEN
            TXT3 = NOMFIC(JJ)(I+3:I+5)
C            WRITE(*,*) ' SBR GETDST. TXT3 = ',TXT3
            READ(TXT3,FMT='(I3)') KDST
            DISTL(J) = DBLE(KDST)/1.D2
            GOTO 2
          ENDIF
          I = I+1
          IF(I.LT.IFI-5) GOTO 1

 2      CONTINUE
C           WRITE(*,*) ' SBR GETDST. I, DIST ',I, DISTL(I)
      ENDDO

      RETURN
      END
