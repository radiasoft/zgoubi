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
      SUBROUTINE HELIXF(X,Y,Z,BR,AK,BO,
     >                                      B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(5,3)
      BN = BO/BR
      AK2 = AK*AK
      AK24 = AK2/4.D0
      AK28 = AK24/2.D0
      XK = X * AK
      CKX = COS(XK)
      SKX = SIN(XK)
      Y2 = Y*Y
      Z2 = Z*Z
      YZ = Y*Z
      Y2Z2 = Y2+Z2
      Y23Z2 = Y2Z2 + 2.D0*Z2
      Z23Y2 = 2.D0*Y2 + Y2Z2

C--------- Components Bx, By, Bz of field
          B(1,1)= -AK*BN*(Y*CKX+Z*SKX)*(1.D0 + AK28*Y2Z2)
          B(1,2)= -BN*((1.D0 + AK28*Z23Y2)*SKX - AK24*YZ*CKX)
          B(1,3)=  BN*((1.D0 + AK28*Y23Z2)*CKX - AK24*YZ*SKX)
 
CC         ... dBx/dX
C          DB(1,1) = 
CC         ... dBx/dY = dBy/dX
C          DB(2,1) = 
CC         ... dBy/dY
C          DB(2,2) = 
CC         .. dBx/dZ = dBz/dX
C          DB(3,1) = 
CC         .. dBy/dZ = dBZ/dY
C          DB(3,2) = 
CC         ... d2Bx/dX2
C          DDB(1,1,1) = 
CC         ... d2Bx/dXdY = d2By/dX2
C          DDB(2,1,1) = 
CC         ... d2Bx/dY2
C          DDB(2,2,1) = 
CC         ... d2By/dY2
C          DDB(2,2,2) = 
CC         .. d2Bx/dXdZ = d2Bz/dX2
C          DDB(3,1,1) = 
CC         .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
C          DDB(3,2,1) = 
CC         .. d2By/dYdZ = d2Bz/dY2
C          DDB(3,2,2) = 
 
      RETURN
      END    
