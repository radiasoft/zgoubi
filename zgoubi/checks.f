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
      SUBROUTINE CHECKS
      use pariz_namelist_interface, only : IZ, ID, MMAP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      INCLUDE 'MXFS.H'

C      IF(MXS.GT.MXD) CALL ENDJOB
C     >('SBR CHECKS ** In MXFS, MXLD, set 10*MXF+MXS-3 .le. MXD',-99)

      IF(IZ.LT.1) 
     >CALL ENDJOB('SBR CHECKS ** In pariz.nml, set IZ to .ge. ',1)

      IF(ID.LT.3) CALL ENDJOB
     >('SBR CHECKS ** In pariz.nml, set ID to .ge. ',3)

      IF(MMAP.LT.1) 
     >CALL ENDJOB(' SBR CHECKS. In pariz.nml, set MMAP to .ge. ',1)

      RETURN
      END
