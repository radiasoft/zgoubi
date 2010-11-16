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
      SUBROUTINE SRDNU(KSC,WF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WF(*)

      PARAMETER(MSAM=1000)
      DIMENSION ANU(MSAM), DNU(MSAM)
      COMMON/FREQ/ ANU, DNU

      MNU = WF(3)
      ANU(1) = WF(1)

      IF(KSC .NE. 3) THEN
C--------- Linear scales

        DNU(1) = ( WF(2) - WF(1) ) / ( MNU -1.D0 )

        DO 1 IN=2,MNU
          ANU(IN) = ANU(IN-1) + DNU(1)
          DNU(IN) =  DNU(1)
 1      CONTINUE
      ELSE
C--------- Log scales

        DNUL = ( DLOG10(WF(2)) - DLOG10(WF(1)) ) / ( WF(3) -1.D0 )
        ANUL = DLOG10(ANU(1))

        DO 2 IN=2,MNU
          ANU(IN) = 10**(ANUL + DNUL)
          ANUL = DLOG10(ANU(IN))
          DNU(IN-1) = ANU(IN) - ANU(IN-1)
C               WRITE(6,*) IN, DNU(IN-1), ANU(IN-1)
 2      CONTINUE
      ENDIF
      RETURN
      END
