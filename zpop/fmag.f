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
      SUBROUTINE FMAG(FNCT,CX,NC,NV,V,
     >                                *)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*), CX(*)

      PARAMETER (MXC=400)
      PARAMETER (MXY=10000)
      DIMENSION X(MXY), Y(MXY)
      EXTERNAL FNCT 

      XL = CX(NC) - CX(1)
      X(1) =  CX(1)
      DX = XL / (MXY - 1.D0)

      DO 3 I = 2, MXY
 3      X(I) = X(I-1) + DX 

      CALL FUNK(FNCT,V,CX(MXC-1),CX(MXC),MXY,NV,X,
     >                                            Y)

      SUM = Y(1)  * DX/2.
      DO 2 I = 2, MXY-1
        DS = Y(I) * DX
        SUM = SUM + DS
 2    CONTINUE
      SUM = SUM + Y(MXY)*DX/2.

      XFB = X(MXY) - SUM

      IF(ABS( CX(MXC-1) - XFB) .GT. DX/2.) THEN
        WRITE(*,FMT='(5X,'' C'',I1,'' =  '',1P,E15.4)')
     >    (I-1,V(I),I=1,NV)
        WRITE(*,*) ' Previous mag face at ',CX(MXC-1)
        WRITE(*,*) ' Adjusted to ', XFB
        WRITE(*,*) ' Busy, matching again...'
        CX(MXC-1) = XFB
        RETURN 1
      ENDIF

      RETURN
      END
