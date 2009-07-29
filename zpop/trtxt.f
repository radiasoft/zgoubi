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
      SUBROUTINE TRTXT(X,Y,TXT,IUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) TXT
      IF(IUN .EQ. 0) THEN
C------- Conversion screen (384/256 Pix)  >  screen (1024/780Pix)
        CALL TRTEXT(SNGL(X*1024./384.), SNGL(Y*780./256.) ,TXT,IUN)
      ELSE
        CALL TRTEXT(SNGL(X),SNGL(Y),TXT,IUN)
      ENDIF
      RETURN
      END
