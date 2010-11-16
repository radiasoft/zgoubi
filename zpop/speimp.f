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
      SUBROUTINE SPEIMP(IUN,YNU,BORNE,U,KT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YNU(*), BORNE(*), U(*)
      INCLUDE 'MAXNTR.H'
      PARAMETER (NTR=NTRMAX*9)
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      IF(YNU(1).GE.BORNE(1) .AND. YNU(1).LE.BORNE(2)) THEN
        XXNU = YNU(1)
      ELSE
        XXNU = 1.D0 - YNU(1)
      ENDIF
      IF(YNU(2).GE.BORNE(3) .AND. YNU(2).LE.BORNE(4)) THEN
        ZZNU = YNU(2)
      ELSE
        ZZNU = 1.D0 - YNU(2)
      ENDIF
      IF(NPTS.EQ.NPTR) 
     >  WRITE(IUN,179) XXNU,ZZNU, (U(I),I=1,3),COOR(KT,6), KT,NPTS
 179    FORMAT(1P,6G14.6,2I6)
      RETURN

      END
