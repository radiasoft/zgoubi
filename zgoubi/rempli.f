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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE REMPLI(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=60) 
      INCLUDE "C.VAR.H"     ! COMMON /VAR/ X(3*MXV),P(MXV)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      DIMENSION F0(6,6)
 
      DO 1 I=1,NV
        K=IR(I)
        L=IS(I)
        IF(K.NE.0) GOTO 10
        CONTINUE
        II1=L/10
        II2=L-10*II1
        IF(M.NE.0) GO TO 9
        CONTINUE
        EXI=F0(1,1)*F0(2,2)-F0(2,1)*F0(2,1)
        EZI=F0(4,4)*F0(5,5)-F0(5,4)*F0(5,4)
        X(I)=F0(II1,II2)

        GOTO 1

 9      CONTINUE
        F0(II1,II2)=X(I)
        IF(II1.GT.2) GO TO 12
        IF(II2.GT.2) GO TO 17
        F0(2,2)=(EXI+F0(2,1)*F0(2,1))/F0(1,1)

        GOTO 17

   12   IF(II1.LT.4.OR.II2.LT.4) GO TO 17
        F0(5,5)=(EZI+F0(5,4)*F0(5,4))/F0(4,4)
   17   CONTINUE

        GOTO 1
 
 10     CONTINUE
        KL=XCOU(I)
        KP=NINT((1D3*XCOU(I)-1D3*KL))
        IF (KL.LT.0) GOTO 2
        IF (KL.GT.0) GOTO 6
        IF(M .EQ. 0)   X(I)=A(K,L)
        IF(M .EQ. 1)      A(K,L)=X(I)
        GOTO 1
 2      CONTINUE
        KL=-KL
        KP=-KP
        CC=A(K,L)+A(KL,KP)
        IF(M.GT.0) GO TO 4
 3      CONTINUE
        X(I)=A(K,L)
        GOTO 1
 4      CONTINUE
        A(K,L) = X(I)
        A(KL,KP) = CC-A(K,L)
        GOTO 1
 6      CONTINUE
        IF(M.LE.0) GO TO 3
        CONTINUE
        A(K,L) = X(I)
        A(KL,KP) = X(I)
 1    CONTINUE
      RETURN
      END
