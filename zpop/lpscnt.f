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
      SUBROUTINE LPSCNT(YM,YPM,U,A,B,XSIGU,NLOG,
     >                                          NCOUNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YM(*), YPM(*), U(*), A(*), B(*), XSIGU(*)
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      NCOUNT = 0
      BETY = B(1)
      ALPY = A(1)
      GAMY = (1.D0+ALPY*ALPY)/BETY
      BETZ = B(2)
      ALPZ = A(2)
      GAMZ = (1.D0+ALPZ*ALPZ)/BETZ
      BETT = B(3)
      ALPT = A(3)
      GAMT = (1.D0+ALPT*ALPT)/BETT
      DO 21 I=1,NPTS
          Y = COOR(I,1)-YM(1)
          YP = COOR(I,2)-YPM(1)
          Z = COOR(I,3)-YM(2)
          ZP = COOR(I,4)-YPM(2)
          T = COOR(I,5)-YM(3)
          TP  = COOR(I,6)-YPM(3)
          IF(GAMY*Y*Y+2.D0*ALPY*Y*YP+BETY*YP*YP.LE.U(1)*XSIGU(1)) THEN 
           IF(GAMZ*Z*Z+2.D0*ALPZ*Z*ZP+BETZ*ZP*ZP.LE.U(2)*XSIGU(2)) THEN 
            IF(GAMT*T*T+2.D0*ALPT*T*TP+BETT*TP*TP.LE.U(3)*XSIGU(3)) THEN
               NCOUNT = NCOUNT + 1
            ENDIF
           ENDIF
          ENDIF
 21   CONTINUE

      WRITE(*,102) U(1),XSIGU(1),U(2),XSIGU(2),U(3),XSIGU(3)
      WRITE(NLOG,102) U(1),XSIGU(1),U(2),XSIGU(2),U(3),XSIGU(3)
 102  FORMAT(/,1P,5X,'3*2-D domains, centered, with emittances :', 
     >/,T25,' Eps_y/pi = ',G12.4,' *',G11.3,' m.rad', 
     >/,T25,' Eps_z/pi = ',G12.4,' *',G11.3,' m.rad',  
     >/,T25,' Eps_t/pi = ',G12.4,' *',G11.3,' s.eV')
      WRITE(*,101) NCOUNT
      WRITE(NLOG,101) NCOUNT
 101  FORMAT(T6,'Number of particles within these 3*2-D domains :',I8,/)

      RETURN
      END
