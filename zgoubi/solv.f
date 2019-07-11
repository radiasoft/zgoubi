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
C  USA
C  -------
C Error, correction B.Bru, 05/2009
C      SUBROUTINE SOLV(A,N,N1,DET,INT)
      SUBROUTINE SOLV(A,N,N1,DET,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RES(20)
      DIMENSION A(20,21),IND(20)
      DET=1.D0
      NN=N1+N
      NT=N
      DO 1 K=1,N
1     IND(K)=K
      DO 2 K=1,N
      RMAX=0.D0
      DO 3 L=K,N
      I1=IND(L)
      R1=ABS(A(I1,K))
      IF(R1.LT.RMAX)GOTO 3
      IMAX=L
      RMAX=R1
3     CONTINUE
      IF(RMAX .EQ. 0.D0)GOTO 100
      I1=IND(K)
      IND(K)=IND(IMAX)
      IND(IMAX)=I1
      INK=IND(K)
      RR=A(INK,K)
      DET=DET*RR
      IF(IMAX .NE. K)DET=-DET
      DO 5 L=K,N
      INL=IND(L)
      IF(INL .EQ. INK)GOTO 5
      RMAX=-A(INL,K)/RR
      DO 6 M=K,NN
6     A(INL,M)=A(INL,M)+RMAX*A(INK,M)
5     CONTINUE
2     CONTINUE
      DET=DET*A(IND(N),N)
      DO 10 L=1,N
      NW=N+1-L
      RES(NW)=A(IND(NW),1+N)/A(IND(NW),NW)
      XX=RES(NW)
      DO 11 M=1,NW
      INM=IND(M)
11    A(INM,N+1)=A(INM,N+1)-A(INM,NW)*XX
10    CONTINUE
      DO 72 K=1,N
72    A(K,N+1)=RES(K)
      RETURN
100   CONTINUE
      DET=0.D0
      RETURN
      END
