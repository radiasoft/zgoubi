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
      SUBROUTINE SERV5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      COMMON//AM(MXX),RM(MXY),ZM(IZ),H(ID,MXX,MXY,IZ),IAMAX,IRMAX,IZMAX
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN
      NX=IAMAX
      NY=IRMAX
      NZ=IZMAX
      DO 1 K=1,NX
1     X(K)=AM(K)
      DO 2 K=1,NY
2     Y(K)=RM(K)
      DO 11 K=1,NZ
11    Z(K)=ZM(K)
C     PRINT 10,(X(K),K=1,NX)
C     PRINT 10,(Y(K),K=1,NY)
C     PRINT 10,(Z(K),K=1,NZ)
10    FORMAT(1X,5E13.6)
177   FORMAT(1X,3HXC=,D20.10,1X,3HYC=,D20.10,3HHX=,D20.10,3HHY=,D20.10)
C     PRINT 10,X(1),X(NX),Y(1),Y(NY)
      RETURN
      END
