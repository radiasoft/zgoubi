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
      SUBROUTINE JCENTR(XT,YT,DX,DY,ICX,ICY,
     >                                      KERK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
C      PARAMETER (MXX=400, MXY=200)
      COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN

      KERK=0

      IF(XT.GE.X(1))GOTO 101

      ICX=2
C      KERK=1
C      PRINT 1000
      GOTO 102

101   CONTINUE
      NNX=NX-1
      DO 1 K=1,NNX
        IF(XT.LT.X(K+1))GOTO 10
1     CONTINUE

      ICX=NX-1
C      KERK=1
C      PRINT 1000
      GOTO 102

10    CONTINUE
      ICX=K
      IF(ABS(XT-X(K)).GT.ABS(XT-X(K+1)))ICX=K+1

102   CONTINUE
      IF(YT.GE.Y(1))GOTO 201

      ICY=2
C      KERK=1
C      PRINT 1000
      GOTO 202

201   CONTINUE
      NNY=NY-1
      DO 2 K=1,NNY
        IF(YT.LT.Y(K+1))GOTO 20
2     CONTINUE
      ICY=NY-1
C      KERK=1
C      PRINT 1000
      GOTO 202

20    CONTINUE
      ICY=K
      IF(ABS(YT-Y(K)).GT.ABS(YT-Y(K+1)))ICY=K+1

202   CONTINUE
      IF(ICX .EQ. 1)ICX=2
      IF(ICY .EQ. 1)ICY=2
      N1X=NX-1
      N1Y=NY-1
      IF(ICX .EQ. NX)ICX=N1X
      IF(ICY .EQ. NY)ICY=N1Y
      DX=XT-X(ICX)
      DY=YT-Y(ICY)
      RETURN
1000  FORMAT(10X,'PARTICLE OUT OF THE MAP ')
      END
