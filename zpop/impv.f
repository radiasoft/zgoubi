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
      SUBROUTINE IMPV(NLOG,NPT,X,Y,YDX,SYDX,IT,KEX,KX,KY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE JUN
      CALL FBGTXT
      WRITE(6 ,FMT='(1P,2G20.12,2G15.5,2I6)') X,Y,YDX,SYDX,NPT,IT
      WRITE(NLOG,FMT='(1P,3E20.12,2I3,A)') X,Y,SYDX,KX,KY,
     >                                        '  X,Y,SYDX,KX,KY'
      WRITE(JUN,FMT='(1P,3E20.12,2I3,A)') X,Y,SYDX,KX,KY,
     >                                        '  X,Y,SYDX,KX,KY'
C      WRITE(NLOG,FMT='(1P,2G20.12,5X,2G18.10,2I6)') X,Y,YDX,SYDX,NPT,IT
      CALL FLUSH2(NLOG,.FALSE.)
      CALL TXTFBG
      RETURN
      ENTRY IMPV2(JUNI)
      JUN = JUNI
      RETURN
      END
