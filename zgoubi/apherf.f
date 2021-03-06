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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      FUNCTION APHERF(XMOY,SIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C     GENERE DES N.A. DANS UNE DISTRIBUTION GAUSSIENNE
C     D'ECART TYPE SIG, CENTREE EN XMOY
C     ------------------------------------------------
      DATA PIS4,R2 / -.7853981634D0, 1.414213562D0/
      IF(RNDM() .GE. 0.5D0) THEN
        SGN=1.D0
      ELSE
        SGN=-1.D0
      ENDIF
      U=1.D0-2.D0*RNDM()
      APHERF= XMOY+SGN*R2*SIG*SQRT(PIS4*LOG(1.D0-U*U))
      RETURN
      END
