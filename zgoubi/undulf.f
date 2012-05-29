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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE UNDULF(BO,X,Z,
     >                           B,DB,DDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      SAVE A, B1
 
      BN=BO*BRI
      P = 1 + (X/A)**8
      P1 = 8.D0 * (X/A)**7 /A
      P2 = 56.D0 * (X/A)**6 /(A*A)
      SBX = SIN(B1*X)
      CBX = COS(B1*X)
      F =  SBX / P  
      F1 =  B1 * CBX / P  -  SBX * P1/(P*P)
      F2 =  - B1*B1 * SBX / P - 2.D0*B1 * CBX * P1/(P*P)  
     >        -  SBX * P2/(P*P) + 2.D0 * SBX * (P1/P)**2/P
      F3 = 0.D0
      F4 = 0.D0
      
      Z2 = Z*Z
      Z3 = Z2*Z
      Z4 = Z3*Z
C  ... Bx, Bz  (By is identically zero)
      B(1,1)= BN * (F1 * Z - F3 * Z3/6.D0)
      B(1,3)= BN * (F - F2 * Z2/2.D0 + F4 * Z4/24.D0)
C  .. dBx/dZ = dBz/dX
      DB(3,1)  = BN * (F1 - F3 * Z2/2.D0)
C  ... dBx/dX
      DB(1,1) = BN * (F2 * Z - F4 * Z3/6.D0)
C  ... d2Bx/dX2
      DDB(1,1,1) = BN * F3 * Z 
C  ... d2Bx/dXdZ = d2Bz/dX2
      DDB(3,1,1) = BN * (F2 - F4 * Z2/2.D0)

      RETURN

      ENTRY UNDUL1(AI,BI)
      A = AI
      B1 = 1.D0/BI
      RETURN
      END
