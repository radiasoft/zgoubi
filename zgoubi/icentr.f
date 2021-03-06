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
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE ICENTR(XT,YT,DX,DY,ICX,ICY,
     >                                      KERK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      INCLUDE "C.SS1.H"     ! COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN

      KERK=0

      IF(XT.GE.X(1))GOTO 101
      ICX=3
C      KERK=1
C      PRINT 1000
      GOTO 102
101   CONTINUE
      NNX=NX-1
      DO 1 K=1,NNX
      IF(XT.LT.X(K+1))GOTO 10
1     CONTINUE
      ICX=NX-2
C      KERK=1
C      PRINT 1000
      GOTO 102
10    CONTINUE
      ICX=K
      IF(ABS(XT-X(K)).GT.ABS(XT-X(K+1)))ICX=K+1
102   CONTINUE
      IF(YT.GE.Y(1))GOTO 201
      ICY=3
C      KERK=1
C      PRINT 1000
      GOTO 202
201   CONTINUE
      NNY=NY-1
      DO 2 K=1,NNY
      IF(YT.LT.Y(K+1))GOTO 20
2     CONTINUE
      ICY=NY-2
      KERK=1
      PRINT 1000
      GOTO 202
20    CONTINUE
      ICY=K
      IF(ABS(YT-Y(K)).GT.ABS(YT-Y(K+1)))ICY=K+1
202   CONTINUE
      IF(ICX .EQ. 1.OR.ICX .EQ. 2)ICX=3
      IF(ICY .EQ. 1.OR.ICY .EQ. 2)ICY=3
      N1X=NX-1
      N2X=NX-2
      N3X=NX-3
      N1Y=NY-1
      N2Y=NY-2
      N3Y=NY-3
      IF(ICX .EQ. NX.OR.ICX .EQ. N1X)ICX=N2X
      IF(ICY .EQ. NY.OR.ICY .EQ. N1Y)ICY=N2Y
      DX=XT-X(ICX)
      DY=YT-Y(ICY)
      RETURN
1000  FORMAT(10X,'PARTICLE OUT OF THE MAP ')
      END
