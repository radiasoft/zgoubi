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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE FLUSH2(IUNIT,BINARY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BINARY, OK
      OK = BINARY   ! just to quiet down the compiler...
C      CHARACTER(80) TXT80
C      IF(IUNIT.EQ.6.OR.IUNIT.EQ.5) RETURN
C      BACKSPACE(IUNIT)
C      IF(BINARY) THEN
C        READ(IUNIT) TXT80
C      ELSE
C        READ(IUNIT,FMT='(A80)') TXT80
C      ENDIF
C 2016, June : I add problem with the above version when adding that in spnprt.
C For this reasoned I chanegd to fortran 2003 style flush :
      IF(IUNIT .GT. 0) FLUSH(IUNIT)
      RETURN
      END
