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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE BXRXYZ(BR,DBR,DDBR,Y,Z,R,ID,B,DB,DDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BR(*),DBR(2,*),DDBR(2,2,*)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
C     -------------------------------------------------------
C     TRANSFORME BR(X,R) EN B(X,Y,Z), AINSI QUE LES DERIVEES.
C     BR = B ou E EN COORD X,R.
C     R=SQRT(Y2+Z2).
C     -------------------------------------------------------
 
      IF(R .NE. 0.D0) THEN
 
        R1=1.D0/R
        DRY=Y*R1
        DRZ=Z*R1
        D2RY=DRZ*DRZ*R1
        D2RYZ=-DRY*DRZ*R1
        D2RZ=DRY*DRY*R1
 
C        Bx, By, Bz
        B(1,1) = BR(1)
        B(1,2) = DRY*BR(2)
        B(1,3) = DRZ*BR(2)
 
        IF(ID .EQ. 0) RETURN
 
C        dBx/dx
        DB(1,1) =  DBR(1,1)
C        dBx/dy=dBy/dx
        DB(2,1) = DRY*DBR(2,1)
C        dBx/dz=dBz/dx
        DB(3,1) = DRZ*DBR(2,1)
C        dBy/dy
        DB(2,2) = R1*( -DRY*B(1,2) + BR(2) + Y*DRY*DBR(2,2) )
C        dBy/dz=dBz/dy
        DB(3,2) = R1*DRZ*( -B(1,2) + Y*DBR(2,2) )
 
        IF(ID .EQ. 1) RETURN
 
C        d2Bx/dx2
        DDB(1,1,1) = DDBR(1,1,1)
C        d2Bx/dxdy=d2By/dx2
        DDB(2,1,1) = DRY*DDBR(2,1,1)
C        d2Bx/dxdz=d2Bz/dx2
        DDB(3,1,1) = DRZ*DDBR(2,1,1)
C        d2Bx/dy2=d2By/dxdy
        DDB(2,2,1) = DRY*DRY*DDBR(2,2,1) + D2RY*DBR(2,1)
C        d2Bx/dydz=d2By/dxdz=d2Bz/dxdy
        DDB(3,2,1) = DRY*DRZ*DDBR(2,2,1) + D2RYZ*DBR(2,1)
C        d2By/dy2
        DDB(2,2,2) = R1*( -D2RY*B(1,2) - 2.D0*DRY*DB(2,2)
     >   + (2.D0*DRY + Y*D2RY)*DBR(2,2) + Y*DRY*DRY*DDBR(2,2,2) )
C        d2By/dydz=d2Bz/dy2
        DDB(3,2,2) = R1*( -D2RYZ*B(1,2) - DRZ*DB(2,2) - DRY*DB(3,2)
     >   + (DRZ + Y*D2RYZ)*DBR(2,2) + Y*DRY*DRZ*DDBR(2,2,2) )
 
      ELSE
 
C        Bx, By, Bz
        B(1,1) = BR(1)
        B(1,2) = 0D0
        B(1,3) = 0D0
        IF(ID .EQ. 0) RETURN
 
C        dBx/dx
        DB(1,1) =  DBR(1,1)
        DB(2,1) = 0D0
        DB(3,1) = 0D0
        DB(2,2) = 0D0
        DB(3,2) = 0D0
        IF(ID .EQ. 1) RETURN
 
        DDB(1,1,1) = DDBR(1,1,1)
        DDB(2,1,1) = 0D0
        DDB(3,1,1) = 0D0
        DDB(2,2,1) = 0D0
        DDB(3,2,1) = 0D0
        DDB(2,2,2) = 0D0
        DDB(3,2,2) = 0D0
 
      ENDIF
 
      RETURN
      END
