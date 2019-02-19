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
      SUBROUTINE COILSF(XX,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM_2.H"     ! COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      PARAMETER (MXCOIL=8)
      INCLUDE "C.COIL.H"     ! COMMON/COIL/ XLI(MXCOIL),RI(MXCOIL),BI(MXCOIL),DIST(MXCOIL),MS

      DO 1 M = 0, MS-1
        M1 = M+1
        IF(M1.GT.MXCOIL)
     >    CALL ENDJOB(' Too  many  coils, max is ',MXCOIL)
        RO = RI(M1)
        BO = BI(M1)
        X=XX-(XS-XE)/2.D0-XE
        CALL SOLEN1(BO,RO)
        CALL SOLENF(X,Y,Z)
 1    CONTINUE

      XL  =-((XS-XE)/2.D0 +X)

      RETURN
      END
