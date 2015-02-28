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
C  -------
      SUBROUTINE eplatf(BN,MPOL,XE,XS,CE,CS,QLE,QLS,QE,QS,X,Y,
     >                                                   E,DE,DDE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MCOEF=6)
      DIMENSION CE(MCOEF), CS(MCOEF)
      DIMENSION QLE(MPOL),QLS(MPOL),QE(MPOL,MCOEF),QS(MPOL,MCOEF)
C      DIMENSION EM(*)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

C-----  dSE/dx=-cosTE=-CTE , dSE/dy = -sinTE =-STE
C-----  dSS/dx=-cosTS=-CTS , dSS/dy = +sinTS =+STS
      SE = ( (XE-X)*CTE - Y*STE )

      IF( SE*SE .LT. 1.D-16) SE=ZERO
 
      IF(SE*SE .LE. XE*XE) THEN
        IF(QLE(1) .NE. ZERO) THEN
          SE=SE/QLE(1)
          CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C     >                                   D7GE,D8GE,D9GE,D10GE)
          DGE  = DGE  * QE(1,1)*BRI 
          D2GE = D2GE * QE(1,2)*BRI 
          D3GE = D3GE * QE(1,3)*BRI 
          D4GE = D4GE * QE(1,4)*BRI 
        ELSE
C--------- QLE=0. => XE=0.
          GE = 1.D0 
          DGE = ZERO
          D2GE= ZERO
          D3GE= ZERO
          D4GE= ZERO
        ENDIF
      ELSE
        IF(SE .GT. ZERO) THEN
          GE= ZERO
        ELSE
          GE = 1.D0
        ENDIF
        DGE = ZERO
        D2GE= ZERO
        D3GE= ZERO
        D4GE= ZERO
      ENDIF
 
      XLS = XLIM-XS
      SS = -( (XS-X)*CTS + Y*STS )

      IF( SS*SS .LT. 1.D-16) SS=ZERO

      IF(SS*SS .LE. XLS*XLS) THEN
        IF(QLS(1) .NE. ZERO) THEN
          SS = SS/QLS(1)
          CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C     >                                   D7GS,D8GS,D9GS,D10GS)
          DGS  = DGS  * QS(1,1)*BRI 
          D2GS = D2GS * QS(1,2)*BRI 
          D3GS = D3GS * QS(1,3)*BRI 
          D4GS = D4GS * QS(1,4)*BRI 
        ELSE
C--------- QLS=0. => XS=0.
          GS = 1.D0
          DGS = ZERO
          D2GS= ZERO
          D3GS= ZERO
          D4GS= ZERO
        ENDIF
      ELSE
        IF(SS .GT. 0.D0) THEN
          GS= ZERO
        ELSE
          GS = 1.D0
        ENDIF
        DGS = ZERO
        D2GS= ZERO
        D3GS= ZERO
        D4GS= ZERO
      ENDIF
 
      CTE2=CTE*CTE
      STE2=STE*STE
      CTS2=CTS*CTS
      STS2=STS*STS

      EY = ( GE+GS-1.D0 )*BN*BRI
 
      EYX=  DGE*CTE - DGS*CTS
      EYY=  DGE*STE + DGS*STS
 
      EYXX= D2GE*CTE2    + D2GS*CTS2
      EYXY= D2GE*CTE*STE - D2GS*CTS*STS
      EYYY= D2GE*STE2    + D2GS*STS2
 
C      EYXXX=  D3GE*CTE2*CTE - D3GS*CTS2*CTS
C      EYXXY=  D3GE*CTE2*STE + D3GS*CTS2*STS
C      EYXYY=  D3GE*CTE*STE2 - D3GS*CTS*STS2
C      EYYYY=  D3GE*STE2*STE + D3GS*STS2*STS
 
c      E(1,1) = 
c      E(1,2) = EY
c      E(1,3) = 0.D0

cC      dEX/dX
c      DE(1,1)= 
cC      dEX/dY, dEY/dX
c      DE(2,1)= 
cC      dEY/dY
c      DE(2,2)= EYY
cC      dEX/dZ, dEZ/dX
c      DE(3,1)=0.D0
cC      dEY/dZ, dEZ/dY
c      DE(3,2)=0.D0

cC      d2EX/dX2 
c      DDE(1,1,1)= 
cC      d2EX/dXdY, d2EY/dX2 
c      DDE(2,1,1)= 
cC      d2EX/dY2, d2EY/dXdY 
c      DDE(2,2,1)= 
cC      d2EY/dY2
c      DDE(2,2,2)= 
cC      d2EX/dXdZ, d2EZ/dX2 
c      DDE(3,1,1)= 
cC      d2EX/dYdZ, d2EY/dXdZ, d2EZ/dXdY
c      DDE(3,2,1) = 
cC      d2EY/dYdZ, d2EZ/dY2 
c      DDE(3,2,2)= 

      RETURN
      END
