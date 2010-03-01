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
      SUBROUTINE SOLBAX(XL,BOSQ,RO2,X,
     >                                BX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BX(*)
      
      X2 = X*X
      U = X2 + RO2
      U2 = U*U
      U3 = U2*U
      SQU = SQRT(U)
      DU = 2.D0 * X
      XL2 = (X-XL)*(X-XL)
      V = XL2 + RO2
      V2 = V*V
      V3 = V2*V
      SQV = SQRT(V)
      DV = 2.D0 * (X-XL)
C------- BX(n)=-d(n-1)B/dX^(n-1)
        BX(1)=BOSQ * ( X/SQU - (X-XL)/SQV )
C           write(88,*) x,bx(1),' solbax'
        BX(2)=BOSQ * ( (U-X2)/(U*SQU) - (V-XL2)/(V*SQV)  )
        BX(3)=BOSQ*( (-4.D0*X*U2 +3.D0*X2*DU*U -U2*DU )/(2.D0*U3*SQU)
     >   - (-4.D0*(XL-X)*V2 +3.D0*XL2*DV*V -V2*DV )/(2.D0*V3*SQV) )
        BX(4)=0.D0
        BX(5)=0.D0
        BX(6)=0.D0
 
      RETURN
      END
