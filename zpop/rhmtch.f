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
C  USA
C  -------
      FUNCTION RHMTCH(NL,MXC,
     >                       NC,CX,CY,SIG,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RHMTCH
      DIMENSION CX(*), CY(*)

      RHMTCH = .FALSE.

      CALL REWIN2(NL,*99)

      NC = 1
 1    CONTINUE
        READ(NL,*,ERR=97,END=98) CX(NC), CY(NC)
        IF(NC .EQ. MXC-2) GOTO 96
        NC = NC + 1
      GOTO 1
 97   WRITE(*,*) ' Stopped  reading. ',NC-1,' data read'

      READ(NL,*,ERR=95,END=95) SIG
      WRITE(*,*) '  I also read Sigma = ', SIG
      READ(NL,*,ERR=95,END=95) NV
      WRITE(*,*) '        degree of the Hermit series = ', NV - 1

      GOTO 95

 96   WRITE(*,*) ' Too many data. Stop reading at ',MXC-2,' data'
      GOTO 95
 98   WRITE(*,*) ' Reached end of file, ',NC-1,' data read'
      GOTO 95

 95   NC = NC -1
      IF(NC .GT. 1) THEN

        RHMTCH = .TRUE.

      ELSE

        WRITE(*,*) ' ***  Not  enough  data  ***' 
        WRITE(*,*) ' -  Check the content of the input data file  -'

      ENDIF

 99   RETURN
      END
