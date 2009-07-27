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
      SUBROUTINE REMPLI(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /INIT/ F0(6,6),F(6,6),EX,EZ,X0,XP0,Z0,ZP0
     >,X1,XP1,Z1,ZP1,PHIX,PHIZ,IF
      PARAMETER (MXV=40) 
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
      DO 1 I=1,NV
        K=IR(I)
        L=IS(I)
        IF(K) 10,7,10
 7      CONTINUE
        II1=L/10
        II2=L-10*II1
        IF(M) 9,8,9
 8      CONTINUE
        EXI=F0(1,1)*F0(2,2)-F0(2,1)*F0(2,1)
        EZI=F0(4,4)*F0(5,5)-F0(5,4)*F0(5,4)
        X(I)=F0(II1,II2)

        GOTO 1

 9      CONTINUE
        F0(II1,II2)=X(I)
        IF(II1-2) 11,11,12
   11   IF(II2-2) 13,13,17
   13   F0(2,2)=(EXI+F0(2,1)*F0(2,1))/F0(1,1)

        GOTO 17

   12   IF(II1-4) 17,18,18
   18   IF(II2-4) 17,19,19
   19   F0(5,5)=(EZI+F0(5,4)*F0(5,4))/F0(4,4)
   17   CONTINUE

        GOTO 1
 
 10     CONTINUE
        KL=NINT(XCOU(I))
        KP=NINT((100.D0*XCOU(I)-100.D0*KL))
        IF(KL) 2,5,6
 5      CONTINUE
        IF(M .EQ. 0)   X(I)=A(K,L)
        IF(M .EQ. 1)      A(K,L)=X(I)
        GOTO 1
C 2      CONTINUE
C        ICOUP=ICOU(I)
C        CC=A(K,1)+A(ICOUP,1)
C        IF(M) 3,3,4
C 3      CONTINUE
C        X(I)=A(K,1)
C        GOTO 1
C 4      CONTINUE
C        A(K,1) = X(I)
C        A(ICOUP,1) = CC-A(K,1)
 2      CONTINUE
        KL=-KL
        KP=-KP
        CC=A(K,L)+A(KL,KP)
        IF(M) 3,3,4
 3      CONTINUE
        X(I)=A(K,L)
        GOTO 1
 4      CONTINUE
        A(K,L) = X(I)
        A(KL,KP) = CC-A(K,L)
        GOTO 1
 6      CONTINUE
        IF(M) 3,3,61
 61     CONTINUE
        A(K,L) = X(I)
        A(KL,KP) = X(I)
 1    CONTINUE
      RETURN
      END
