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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE XROTB(A,B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C----- composantes du champ dans les axes tournes de A (sens trigo)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)

      DYY=COS(A)
      DZZ=DYY
      DYZ=SIN(A)
      DZY=-DYZ
      DYY2=DYY*DYY
      DZZ2=DZZ*DZZ
      DYZ2=DYZ*DYZ

      B12=B(1,2)*DYY+B(1,3)*DZY
      B(1,3)=B(1,2)*DYZ+B(1,3)*DZZ
      B(1,2)=B12

      DB33  = - ( DB(1,1) + DB(2,2) )

C        dBx/dy=dBy/dx
        DB21 = DB(2,1)*DYY+DB(3,1)*DZY
C        dBx/dz=dBz/dx
        DB(3,1) = DB(2,1)*DYZ+DB(3,1)*DZZ
C        dBy/dy
        DB22 = DB(2,2)*DYY2+2.D0*DB(3,2)*DYY*DZY+DB33*DYZ2
C        dBy/dz=dBz/dy
        DB32 = DB(2,2)*DYZ*DYY+DB(3,2)*(DZZ*DYY+DYZ*DZY)+(DB33*DZZ*DZY)

        DB(2,1)=DB21
        DB(2,2)=DB22
        DB(3,2)=DB32

      DDB331 = -(DDB(1,1,1) + DDB(2,2,1))
      DDB332 = -(DDB(2,1,1) + DDB(2,2,2))
      DDB333 = -(DDB(3,1,1) + DDB(3,2,2))

C        d2Bx/dxdy=d2By/dx2
        DDB211 = DDB(2,1,1)*DYY+DDB(3,1,1)*DZY
C        d2Bx/dxdz=d2Bz/dx2
        DDB311 = DDB(2,1,1)*DYZ+DDB(3,1,1)*DZZ
C        d2Bx/dy2=d2By/dxdy
        DDB221 = DDB(2,2,1)*DYY2+2.D0*DDB(3,2,1)*DZY*DYY
     >   +DDB331*DYZ2
C        d2Bx/dydz=d2By/dxdz=d2Bz/dxdy
        DDB321 = DDB(2,2,1)*DYY*DYZ+DDB(3,2,1)*(DZY*DYZ+DYY*DZZ)
     >   +DDB331*DZY*DZZ
C        d2By/dy2
        DDB222 = ( DDB(2,2,2)*DYY+3.D0*DDB(3,2,2)*DZY )*DYY2
     >   + ( 3.D0*DDB332*DYY+DDB333*DYZ )*DYZ2
C        d2By/dydz=d2Bz/dy2
        DDB322=DDB(2,2,2)*DYZ*DYY2+DDB(3,2,2)*DYY*(2.D0*DZY*DYZ+DYY*DZZ)
     >   + DDB332*DZY*(DZY*DYZ+2.D0*DYY*DZZ)+DDB333*DYZ2*DZZ

        DDB(2,1,1) =DDB211
        DDB(3,1,1) =DDB311
        DDB(2,2,1) =DDB221
        DDB(3,2,1) =DDB321
        DDB(2,2,2) =DDB222
        DDB(3,2,2) =DDB322

      DO 1 M=1,3
        DO 1 L=1,3
          DO 1 J=1,3
            D3BX(J,L,M)=0.D0
            D3BY(J,L,M)=0.D0
            D3BZ(J,L,M)=0.D0
            DO 1 I=1,3
              D4BX(I,J,L,M)=0.D0
              D4BY(I,J,L,M)=0.D0
              D4BZ(I,J,L,M)=0.D0
1     CONTINUE

      RETURN
      END
