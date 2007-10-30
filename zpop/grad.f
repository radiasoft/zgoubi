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
      SUBROUTINE GRAD(FNCT,DFNCT,C,CX,CY,NC,NV,FCT,RO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*),RO(*),CX(*),CY(*),FCT(*)
      PARAMETER(MXV=10)
      PARAMETER(MXC=400)
      EXTERNAL FNCT,DFNCT

      XFB= CX(MXC-1)
      XLAMB= CX(MXC)

      DO 10 J=1,NV      
10      RO(J)=0.D0

      DO 20 I = 1, NC
        T = ( CX(I) - XFB ) / XLAMB
        DO 20 J = 1, NV
          DS = (FCT(I)-CY(I)) * DFNCT(J,C,T)
20        RO(J) = RO(J) + DS

      DO 21 J = 1, NV
        RO(J) = 2.D0 * RO(J)
 21   CONTINUE

      RETURN
      END               
