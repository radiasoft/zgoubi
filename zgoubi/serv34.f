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
      SUBROUTINE SERV34
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARIZ.H'
      COMMON/SS1/X(MXX),Y(MXY),Z(IZ),JY(25),JX(25),JZ(25),NX,NY,NZ,NN
      COMMON/SSS/EE(15,25)
      DIMENSION A(20,21),IND(20),D(25)
      DO 100 KK=1,NN
        DO 101 LL=1,NN
101     D(LL)=0.D0
      D(KK)=1.D0
      IXC=2
      IYC=2
      DO 1 K=1,6
        DO 1 L=1,6
        CC=0D0
        DO 2 M=1,NN
          IX=IXC+JX(M)
          IY=IYC+JY(M)
          XX=X(IX)-X(IXC)
          YY=Y(IY)-Y(IYC)
 2        CC=CC+GG(K,XX,YY)*GG(L,XX,YY)
        A(K,L)=CC
 1      CONTINUE
      DO 3 K=1,6
        CC=0.D0
        DO 4 M=1,NN
          IX=IXC+JX(M)
          IY=IYC+JY(M)
          XX=X(IX)-X(IXC)
          YY=Y(IY)-Y(IYC)
4         CC=CC+GG(K,XX,YY)*D(M)
        A(K,7)=CC
3       CONTINUE
      CALL SOLV(A,6,1,DET,IND)
      EE(1,KK)=A(1,7)
      EE(2,KK)=A(2,7)
      EE(3,KK)=A(3,7)
      EE(4,KK)=A(4,7)
      EE(5,KK)=A(5,7)
      EE(6,KK)=A(6,7)
100   CONTINUE
      RETURN
      END
