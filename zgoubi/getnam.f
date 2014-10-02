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
      SUBROUTINE GETNAM(LUN,MXDST,
     >                            NOMFIC,NBDST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) NOMFIC(*)
      INTEGER DEBSTR,FINSTR
      I = 1
 1    CONTINUE
        II = I+1
        READ(LUN,*,END=98,ERR=99) NOMFIC(I)
        NOMFIC(I) = NOMFIC(I)(DEBSTR(NOMFIC(I)):FINSTR(NOMFIC(I)))
        READ(LUN,*,END=98,ERR=99) NOMFIC(II)
        NOMFIC(II) = NOMFIC(II)(DEBSTR(NOMFIC(II)):FINSTR(NOMFIC(II)))
        IF(I.EQ.MXDST) GOTO 9
        I = I+2
        GOTO 1

 9    CONTINUE
      WRITE(6,*) ' SBR GETNAM. Reached max number of map files allowed'
      NBDST = I/2
      RETURN

 98   CONTINUE
      NBDST = I/2
      RETURN

 99   STOP ' SBR GETNAM. Error while reading '
      END
