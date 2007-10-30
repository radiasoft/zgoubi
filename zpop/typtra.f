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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE TYPTRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------
C     LINE TYPE
C     ---------

      CALL LINTYR(LTYP)

 21   CONTINUE
      CALL HOMCLR
      
      WRITE(6,100) LTYP
 100  FORMAT(//,3X,60(1H*),//,20X,' LINE  TYPE :',//
     X,20X,' ( Now =  ',I2,' )',//
     1 ,9X,'  1     solid         ',/ 
     2 ,9X,'  2    -.-.-.-        ',/ 
     3 ,9X,'  3    .......        ',/ 
     4 ,9X,'  4    .......        ',/ 
     5 ,9X,'  5    -------        ',/ 
     6 ,9X,'  6    -.-.-.-        ',/  
     7 ,9X,'  7    -------        ',/ 
     8 ,9X,'  8    -.-.-.-        ',/  
     9 ,9X,'  9       .           ',/  
     X ,9X,' 10       +           ',/  
     1 ,9X,' 11       *           ',/  
     2 ,9X,' 12       o           ',/  
     3 ,9X,' 13       x           ',///)

      WRITE(6,102) ' * Option  number : '
 102  FORMAT(A20,$)
      READ(5,FMT='(I2)',ERR=21) LTYP 
      IF( LTYP .LT. 1  .OR.  LTYP .GT. 13 ) GOTO 21
      CALL LINTYW(LTYP)
      RETURN
      END
