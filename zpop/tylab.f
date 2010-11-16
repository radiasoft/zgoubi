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
      FUNCTION TYLAB(KX,KY)
C--------------------------------------------------
C     TYLAB is .TRUE. if variables to be plotted are
C        both of type  Xlab, Ylab or Zlab
C--------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TYLAB,TYLABR,LAB
      SAVE LAB
c      TYLAB = (KX .EQ. 42 .OR. KX .EQ. 44 .OR. KX .EQ. 48 
c     >    .OR. KX .EQ. 62 .OR. KX .EQ. 68)
c     >       .AND. (KY .EQ. 42 .OR. KY .EQ. 44 .OR. KY .EQ. 48
c     >    .OR. KY .EQ. 62 .OR. KY .EQ. 68)
      TYLAB = (KX .EQ. 42 .OR. KX .EQ. 44 .OR. KX .EQ. 48)
     >  .AND. (KY .EQ. 42 .OR. KY .EQ. 44 .OR. KY .EQ. 48)
      LAB=TYLAB
      IF(LAB) WRITE(6,FMT='(/,A)') ' Lab. coordinates. '
      RETURN
      ENTRY TYLABR()
      TYLABR = LAB
      RETURN
      END
