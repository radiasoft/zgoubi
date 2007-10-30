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
      SUBROUTINE MNFILE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

 21   CONTINUE
      CALL HOMCLR

      WRITE(6,100) 
 100  FORMAT(//,3X,60(1H*),//,20X,' MENU :',//
     1 ,9X,'  1    FICHIERS DATA         ( .dat )   ',/
     2 ,9X,'  2    FICHIERS RESULTAT     ( .res )   ',/
     3 ,9X,'  3    FICHIERS TRAJECTOIRES ( .plt )   ',/
     3 ,9X,'          (OPTION  ''IL=2'')     ',/
     4 ,9X,'  4    FICHIERS FAISCEAU     ( .fai )   ',/
     4 ,9X,'          (MOT-CLE  ''FAISCNL'')  ',/
     5 ,9X,'  5    FICHIERS SPIN         ( .spn )   ',/
     9 ,9X,'  9        EXIT THIS MENU        ',///,3X,60(1H*),/)

      WRITE(6,105) ' * Option  number : '
 105  FORMAT(A20,$)
      READ(5,101,ERR=21) IOPT
 101  FORMAT(I2)
      GOTO ( 1, 2, 3, 4, 5,21,21,21,99) IOPT
      GOTO 21

 1    CONTINUE                                          
         CALL MNFICH('.dat') 
         GOTO 21

 2    CONTINUE                                          
         CALL MNFICH('.res') 
         GOTO 21
 
 3     CONTINUE
         CALL MNFICH('.plt') 
         GOTO 21       
 
 4     CONTINUE
         CALL MNFICH('.fai') 
         GOTO 21       

 5     CONTINUE
         CALL MNFICH('.spn') 
         GOTO 21       

 99    RETURN     
       END
