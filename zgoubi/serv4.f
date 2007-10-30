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
      SUBROUTINE SERV4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'PARIZ.H'
C      PARAMETER (MXX=400, MXY=200)
      COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN
      COMMON/SSS/EE(15,25)
      DIMENSION A(20,21),IND(20),D(25)

      DO 100 KK=1,NN
      DO 101 LL=1,NN
101   D(LL)=0.D0
      D(KK)=1.D0
      IXC=3
      IYC=3
      DO 1 K=1,15
      DO 1 L=1,15
      CC=0.D0
      DO 2 M=1,NN
      IX=IXC+JX(M)
      IY=IYC+JY(M)
      XX=X(IX)-X(IXC)
      YY=Y(IY)-Y(IYC)
2     CC=CC+FGF(K,XX,YY)*FGF(L,XX,YY)
      A(K,L)=CC
1     CONTINUE
      DO 3 K=1,15
      CC=0.D0
      DO 4 M=1,NN
      IX=IXC+JX(M)
      IY=IYC+JY(M)
      XX=X(IX)-X(IXC)
      YY=Y(IY)-Y(IYC)
4     CC=CC+FGF(K,XX,YY)*D(M)
      A(K,16)=CC
3     CONTINUE
      CALL SOLV(A,15,1,DET,IND)
      EE(1,KK)=A(1,16)
      EE(2,KK)=A(2,16)
      EE(3,KK)=A(3,16)
      EE(4,KK)=A(4,16)
      EE(5,KK)=A(5,16)
      EE(6,KK)=A(6,16)
      EE(7,KK)=A(7,16)
      EE(8,KK)=A(8,16)
      EE(9,KK)=A(9,16)
      EE(10,KK)=A(10,16)
      EE(11,KK)=A(11,16)
      EE(12,KK)=A(12,16)
      EE(13,KK)=A(13,16)
      EE(14,KK)=A(14,16)
      EE(15,KK)=A(15,16)
100   CONTINUE
      RETURN
      END
