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
       SUBROUTINE BAXBXR(BX,R,R2,B,DB,DDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BX(*),B(*),DB(2,*),DDB(2,2,*)
C     -------------------------------------------------
C     Taylor series yielding  B(X,R) & derivatives as functions 
C     of R, computed from  dBX/dX at R=0 (= on axis).
C     B = B or E in X,R  coordinates. 
C     BX = field on X axis = BX(X,R=0).
C     BX(n+1)=dnBX/dXn.
C     -------------------------------------------------
 
C----- BX(X,R), BR(X,R)
      B(1)= BX(1)+( -.25D0*BX(3) + .015625D0*BX(5)*R2 )*R2
      B(2)= (-0.5D0*BX(2)+(.0625D0*BX(4)-R2*BX(6)/384.D0)*R2 )*R
 
C----- DERIVEES 1-ERES DE B
C     INDICE 1 = ELEMENT DIFFERENTIEL: dX, dR
C     INDICE 2 DU TABLEAU = COMPOSANTE: BX, BR
 
C      dBX/dX
      DB(1,1) = BX(2)+(-.25D0*BX(4) + .015625D0*BX(6)*R2)*R2
C      dBX/dR=dBR/dX
      DB(2,1) = ( -.5D0*BX(3) + .0625D0*R2*BX(5) )*R
C      dBR/dR
      DB(2,2) = (-.5D0*BX(2)+ ( .1875D0*BX(4) -5.D0*R2*BX(6)/384D0 )*R2)
 
C----- DERIVEES 2-EMES DE B
C     INDICES 1,2 = ELEMENT DIFFERENTIEL: dX, dR
C     INDICE 3 DU TABLEAU = COMPOSANTE: BX, BR
 
C      d2BX/dX2
      DDB(1,1,1) = BX(3) - .25D0*R2*BX(5)
C      d2BX/dXdR=d2BX/dRdX=d2BR/dX2
      DDB(2,1,1) = ( -.5D0*BX(4) + .0625D0*R2*BX(6) )*R
C      d2BR/dXdR=d2BX/dR2
      DDB(2,2,1) = -.5D0*BX(3) + .1875D0*R2*BX(5)
C      d2BR/dR2
      DDB(2,2,2) = (.375D0*BX(4) -5.D0*R2*BX(6)/96D0 )*R
 
      RETURN
      END
