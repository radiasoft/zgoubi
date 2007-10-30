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
      SUBROUTINE UNDULFA(BN,MPOL,DLE,DLS,DE,DS,RTB,X,Y,
     >                                                BZ0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DLE(MPOL),DLS(MPOL),DE(MPOL,10),DS(MPOL,10),RTB(MPOL)
      DIMENSION BZ0(5,5)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/RIGID/ BORO,DPREF,DP,BR
 
      SE = XE-X
      IF( SE*SE .LT. 1.D-8) SE=ZERO
       IF(SE*SE .LE. XE*XE) THEN
        IF(DLE(1) .NE. ZERO) THEN
          SE=SE/DLE(1)
          CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C     >                                   D7GE,D8GE,D9GE,D10GE)
          DGE  = DGE  * DE(1,1)/BR
          D2GE = D2GE * DE(1,2)/BR
          D3GE = D3GE * DE(1,3)/BR
          D4GE = D4GE * DE(1,4)/BR
        ELSE
C--------- DLE=0. => XE=0.
          GE = 1D0 
          DGE = ZERO
          D2GE= ZERO
          D3GE= ZERO
          D4GE= ZERO
        ENDIF
      ELSE
        IF(SE .GT. ZERO) THEN
          GE= ZERO
        ELSE
          GE = 1D0
        ENDIF
        DGE = ZERO
        D2GE= ZERO
        D3GE= ZERO
        D4GE= ZERO
      ENDIF
 
      XLS = XLIM-XS
      SS = -(XS-X)
      IF( SS*SS .LT. 1.D-8) SS=ZERO 
      IF(SS*SS .LE. XLS*XLS) THEN
        IF(DLS(1) .NE. ZERO) THEN
          SS = SS/DLS(1)
          CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C     >                                   D7GS,D8GS,D9GS,D10GS)
          DGS  = DGS  * DS(1,1)/BR
          D2GS = D2GS * DS(1,2)/BR
          D3GS = D3GS * DS(1,3)/BR
          D4GS = D4GS * DS(1,4)/BR
        ELSE
C--------- DLS=0. => XS=0.
          GS = 1D0
          DGS = ZERO
          D2GS= ZERO
          D3GS= ZERO
          D4GS= ZERO
        ENDIF
      ELSE
        IF(SS .GT. 0.D0) THEN
          GS= ZERO
        ELSE
          GS = 1D0
        ENDIF
        DGS = ZERO
        D2GS= ZERO
        D3GS= ZERO
        D4GS= ZERO
      ENDIF
 
      FAC = BN/BR
      UK2 = UK*UK 
      BZ = FAC * SIN(UK * (X-XE))
      BZX=  FAC * UK * COS(UK * (X-XE))
      BZXX= - UK2 * BZ
      BZXXX= - UK2 * BZX
 
C Due to reminicences from ancient times
      BZ0(1,1)=BZ
      BZ0(2,1)=BZX
      BZ0(3,1)=BZXX
      BZ0(4,1)=BZXXX

      RETURN
      ENTRY UNDUL2(XL,NP)
      UK = 2.D0*PI/(XL/NP)
      RETURN
      END
