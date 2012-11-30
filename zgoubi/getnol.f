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
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE GETNOL(KLEY,
     >                       NUML)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (KSIZ=10)
      CHARACTER*(KSIZ)  KLEY, KLE
      KLE = '*******'
      NUML = 0
      DO WHILE ( KLEY .NE. KLE .AND. NUML.LE.MXL )
        NUML = NUML+1
        CALL ZGKLE(IQ(NUML)
     >                     ,KLE)
      ENDDO
      IF(.NOT. (KLEY .EQ. KLE)) THEN
        CALL ENDJOB('SBR getnol : no such key #',IQ(NUML))
      ENDIF       
      RETURN
      END
