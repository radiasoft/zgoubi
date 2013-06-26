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
      SUBROUTINE MNZGRA(IOPT,wrkDir)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character(*) wrkDir
      integer debstr, finstr

c      write(*,*) ' Working directory : ', 
c     >   wrkDir(debstr(wrkDir):finstr(wrkDir))
c      read(*,*)

 21   CONTINUE
      CALL HOMCLR

      WRITE(6,100) '[...]'//wrkDir(finstr(wrkDir)-60:finstr(wrkDir))
 100  FORMAT(//,3X,60('*'),//,20X,' MENU :', //
     > ,4x,' Working directory now : ',/,5x,A,  //
c     1 ,9X,'  1    Run   Zgoubi   ',/
c     2 ,9X,'  2    Batch Zgoubi   ',/
c     3 ,9X,'  3    List files ( .dat, .res, .plt, .fai, .spn ) ',/
     7 ,9X,'  7    Menu  for  direct  plots  from',/
     7 ,9X,'          .plt, .fai, .spn  type  zgoubi output files ',/
     > ,9X,'  ',/
     8 ,9X,'  8    Menu  for other  Analysis/Plot  on ',/
     7 ,9X,'          .plt, .fai, .spn  type  zgoubi output files ',/
     > ,9X,'  ',/
     X ,9X,' 10    Exit  zpop  ',/
c     > ,9X,'  ',/
c     > ,9X,'  --------',/
c     > ,9X,'  AGS MODEL SPECIFIC : ',/
c     > ,9X,'  ',/
c     1 ,9X,' 11    ZgoubiFromSnaprampCmd ',/
c     1 ,9X,' 12    Tune scan from sdds files '
     >,///,3X,60('*'),/)

      WRITE(6,105) ' * Option  number : '
 105  FORMAT(A20)
      READ(5,101,ERR=21) IOPT
 101  FORMAT(I2)
      IF(IOPT.GE. 7 .OR. IOPT.LE. 10) THEN
        RETURN 
      ELSE
        GOTO 21
      ENDIF
      END
