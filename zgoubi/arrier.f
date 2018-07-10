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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE ARRIER(S)

      CHARACTER(*) S

      INTEGER MODE

      CHARACTER(1) ANSI_CR
      CHARACTER(2) ANSI_CSI
      CHARACTER(3) ANSI_CUU_1
      INTEGER TTYMODE
      SAVE TTYMODE,ANSI_CR,ANSI_CUU_1

      IF (TTYMODE.EQ.1) THEN
         WRITE(6,'(''+'',A,A)') ANSI_CUU_1//ANSI_CR,S
      ELSE
         WRITE(6,'(''+'',A)') S
      ENDIF
      
      RETURN

      ENTRY INIARR(MODE)

      TTYMODE=MODE
      ANSI_CR    = CHAR(13)
      ANSI_CSI   = CHAR(27)//'['
      ANSI_CUU_1 = ANSI_CSI//'A'

      RETURN

      END
