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
      SUBROUTINE SERV5
      use xyzhc_interface, only : XH, YH, ZH, IXMA, JYMA, KZMA
      use c_ss1_interface, only : X, Y, Z, NX, NY, NZ
      USE dynhc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CHAVE.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      NX=IXMA
      NY=JYMA
      NZ=KZMA
      DO 1 K=1,NX
1     X(K)=XH(K)
      DO 2 K=1,NY
2     Y(K)=YH(K)
      DO 11 K=1,NZ
11    Z(K)=ZH(K)
C     PRINT 10,(X(K),K=1,NX)
C     PRINT 10,(Y(K),K=1,NY)
C     PRINT 10,(Z(K),K=1,NZ)
C     PRINT 10,X(1),X(NX),Y(1),Y(NY)
      RETURN
      END
