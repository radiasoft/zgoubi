C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SYNBOX(N,TTA1,S1,Y1,TTA,S,Y,AL,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARCST.H'

      TTA=TTA1+PISUR2
      S=S1+W*COS(TTA)
      Y=Y1+W*SIN(TTA)
      CALL CALTRI(N,S,Y)

      TTA=TTA1
      S=S+AL*COS(TTA)
      Y=Y+AL*SIN(TTA)
      CALL CALTRI(N,S,Y)

      TTA=TTA1-PISUR2
      S=S+2.*W*COS(TTA)
      Y=Y+2.*W*SIN(TTA)
      CALL CALTRI(N,S,Y)

      TTA=TTA1-PI
      S=S+AL*COS(TTA)
      Y=Y+AL*SIN(TTA)
      CALL CALTRI(N,S,Y)
 
      TTA=TTA1-3.*PISUR2
      S=S+W*COS(TTA)
      Y=Y+W*SIN(TTA)
      CALL CALTRI(N,S,Y)

      CALL SYNAXE(N,TTA1,S,Y,TTA,S,Y,AL)

      RETURN
      END
