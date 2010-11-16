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
      SUBROUTINE SRMODL(NLOG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GPH(10), GPS(10)
      CHARACTER REP*1
      DATA KC / 1 /
      DATA KP, NPH, NPS / 1, 1, 1 /
      DATA GPH / 100.D0, 100.D0, 4.D0, 7 * 0.D0 /
      DATA GPS / 10 * 0.D0 / 

      GOTO 21
      
 20   CONTINUE      
      CALL FBGTXT
      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
      READ(5,FMT='(A1)',ERR=20) REP 

 21   CONTINUE
      CALL FBGTXT
C      CALL HOMCLR

      WRITE(6,104) 
 104  FORMAT(//,3X,60('-')
     >,//,10X,' Synchrotron Radiation, E-field Models',/
     1,/,5X,' 1  Chose  sigma  or  pi  component '
     2,/,5X,' 2  Give  Gamma.Phi  and  Gamma.Psi '
     3,/,5X,' 3  Give  physical  parameters'
     7,/,5X,' 7  ** Plot **'
     9,/,5X,' 9  Exit  this  menu '
     >,/,3X,60('-'),//)

      WRITE(6,100)
 100  FORMAT('$  YOUR  CHOICE : ')
      READ(5,108,ERR=21) IOPT
108   FORMAT(I2)
      GOTO ( 1, 2, 3,21,21,21, 7,21,99) IOPT  
      GOTO 21
 
 1    CONTINUE
      GOTO 21

 2    CONTINUE
 24   WRITE(6,*) ' How many curves on the same plot:'
      READ(5,*,err=24) ncu
      WRITE(6,*) ' parameter = Gamma*Phi (1)  or Gamma*Psi (2):'
      READ(5,*) kp
      if(kp .eq. 1) then
        nph = ncu
        nps = 1
        WRITE(6,*) ' Give the ',nph,' values for Gamma*Phi:'
        READ(5,*) (gph(i), i=1, nph)
        WRITE(6,*) ' Give the value for Gamma*Psi:'
        READ(5,*) gps(1)
      else
        nph = 1
        nps = ncu
        WRITE(6,*) ' Give the ',nps,' values for Gamma*Psi:'
        do 23 i=1, nps
 23        READ(5,*) gps(i)
        WRITE(6,*) ' Give the value for Gamma*Phi:'
        READ(5,*) gph(1)
      endif
      GOTO 21

 3    CONTINUE
      GOTO 21

 7    CONTINUE
        CALL SRMDL1(NLOG,KC,NPH,GPH,NPS,GPS)
      GOTO 20
      
 99   RETURN
      END
