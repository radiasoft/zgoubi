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
      SUBROUTINE MNZGRA(IOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

 21   CONTINUE
      CALL HOMCLR

      WRITE(6,100) 
 100  FORMAT(//,3X,60('*'),//,20X,' MENU :',///
C     1 ,9X,'  1    RUN   ZGOUBI   ',/
C     2 ,9X,'  2    BATCH ZGOUBI   ',/
C     3 ,9X,'  3    FILES ( .dat, .res, .plt, .fai, .spn ) ',/
     7 ,9X,'  7    Menu  for  direct  plots  from',/
     7 ,9X,'          .plt, .fai, .spn  type  zgoubi output files ',/
     > ,9X,'  ',/
     8 ,9X,'  8    Menu  for other  Analysis/Plot  ',/
     > ,9X,'  ',/
     X ,9X,' 10    Exit  zpop  ',/
     > ,9X,'  ',/
C     1 ,9X,' 11    SHOW QUEUE '
     >,///,3X,60('*'),/)

      WRITE(6,105) ' * Option  number : '
 105  FORMAT(A20)
      READ(5,101,ERR=21) IOPT
 101  FORMAT(I2)
      IF(IOPT.GE. 1 .OR. IOPT.LE. 11) THEN
        RETURN 
      ELSE
        GOTO 21
      ENDIF
      END
