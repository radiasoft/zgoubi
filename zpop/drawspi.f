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
      subroutine drawspi(tta0,r0,rp1,rp2,xi,lunOut)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=100)
      pi = 4.d0 * atan(1.d0)
      b = 1.d0 / tan(xi*pi/180.d0)

C Limit angles corresponding to limit radii rp1, rp2
        tta1 =  log(rp1/r0) / b
        tta2 =  log(rp2/r0) / b

Compute spiral EFB
        DO 100 I=0,N
           TTA = ((N-I)*tta1 + I*tta2)/N
           r = r0 *exp(b * tta)
           x = r *cos(tta + tta0)  
           y = r *sin(tta + tta0)
           curv = exp( - b * tta)/ (r0 * sqrt(1 + b*b))
           write(lunOut,*) x, y, tta, curv
 100    CONTINUE
      
      return
      end
