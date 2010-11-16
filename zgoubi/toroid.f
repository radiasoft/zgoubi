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
      SUBROUTINE TOROID(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CHAMP/ BZ0(5,5), EZ0(5,5)
      COMMON/CHAVE/ B(15),V(5,3),E(5,3)
      COMMON/DDBXYZ/ DB(9),DDB(27)
      COMMON/D3BXYZ/ D3BX(27), D3BY(27), D3BZ(27)
      COMMON/D4BXYZ/ D4BX(81) ,D4BY(81) ,D4BZ(81)
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
      COMMON/MARK/ KART,KALC,KERC,KUASEX
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      EQUIVALENCE (R1,EB1),(R2,EB2),(R3,EG1),(R4,EG2)
 
      BN=BO*BRI
          XX = X
          X = ABS(X)
 
          IF(Y*Y .LE. RO*RO) THEN
 
            IF    (EN .EQ. 1.D0) THEN
              IF(X .LT. R1) THEN
                BZ    = 0D0
                BZX   = 0D0
                BZXX  = 0D0
              ELSEIF(X .LT. R2) THEN
                R12=R1*R1
C               ... BN*R2 = B(R2)*R2 = Mu0*N*I/(2*PI*ALPHA)
                FAC = BN*R2/(R2*R2-R12)
                BZ    = FAC * ( X - R12/X)
C                BZX   = FAC * ( 1.D0- R12/X/X )
                BZX   = FAC * ( 1.D0+ R12/X/X )
                BZXX  = -2.D0 * FAC * R12/X/X/X
              ELSEIF(X .LT. R3) THEN
                FAC = BN*R2
                BZ    = FAC / X
                BZX   = -BZ/X
                BZXX  = -2.D0*BZX/X
              ELSEIF(X .LT. R4) THEN
                R42=R4*R4
                FAC = -BN*R2/(R42-R3*R3)
                BZ    = FAC * ( X - R42/X )
C                BZX   = FAC * ( 1.D0- R42/X/X )
                BZX   = FAC * ( 1.D0+ R42/X/X )
                BZXX  = -2.D0 * FAC * R42/X/X/X
              ELSEIF(X .GE. R4) THEN
                BZ    = 0D0
                BZX   = 0D0
                BZXX  = 0D0
              ENDIF
 
            ELSEIF(EN .EQ. 2.D0) THEN
              IF(X .LT. R1) THEN
                BZ    = 0D0
                BZX   = 0D0
              ELSEIF(X .LT. R2) THEN
                BZX   = BN / (R2-R1)
                BZ    = BZX* ( X  - R1 )
              ELSEIF(X .LT. R3) THEN
                BZ    = BN
                BZX   = 0D0
              ELSEIF(X .LT. R4) THEN
                BZX = -BN/(R4-R3)
                BZ   = BZX * (X - R4 )
              ELSEIF(X .GE. R4) THEN
                BZ    = 0D0
                BZX   = 0D0
              ENDIF
            ENDIF
C           ... ENDIF EN=1,2
 
            X=XX
            IF(X.LT. 0.D0) THEN
              BZ = -BZ
              BZX = -BZX
              BZXX = -BZXX
            ENDIF
 
          ELSE
            BZ    = 0D0
            BZX   = 0D0
            BZXX  = 0D0
          ENDIF
C         ... ENDIF Y**2 <> RO**2
      RETURN
      END
