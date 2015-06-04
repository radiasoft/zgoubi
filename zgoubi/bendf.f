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
      SUBROUTINE BENDF(BN,MPOL,XE,XS,CE,CS,DLE,DLS,DE,DS,RTB,X,Y,
     >                                                           BZ0)
C see FM CeeRainer 11/2003. QCE in CHAMC
C      SUBROUTINE BENDF(BN,MPOL,DLE,DLS,DE,DS,RTB,X,Y,
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MCOEF=6)
      DIMENSION CE(MCOEF), CS(MCOEF)
      DIMENSION DLE(MPOL),DLS(MPOL),
     >  DE(MPOL,MCOEF),DS(MPOL,MCOEF),RTB(MPOL)
      DIMENSION BZ0(5,*)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

      SAVE IOPIN, IOPOU, DISTE, DISTS
      SAVE G1, G2

      DATA G1, G2 /0d0,0d0/
             
      CTE=RTB(1)
      STE=RTB(2)
      CTS=RTB(4)
      STS=RTB(5)
C-----  dSE/dx=-cosTE=-CTE , dSE/dy = -sinTE =-STE
C-----  dSS/dx=-cosTS=-CTS , dSS/dy = +sinTS =+STS
      SE = ( (XE-X)*CTE - Y*STE )

      IF( SE*SE .LT. 1.D-16) SE=ZERO
 
      IF(SE*SE .LE. XE*XE) THEN
        IF(DLE(1) .NE. ZERO) THEN
          SE=SE/DLE(1)
          CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C     >                                   D7GE,D8GE,D9GE,D10GE)
          DGE  = DGE  * DE(1,1)*BRI 
          D2GE = D2GE * DE(1,2)*BRI 
          D3GE = D3GE * DE(1,3)*BRI 
          D4GE = D4GE * DE(1,4)*BRI 
          IF(IOPIN.LT.0) THEN
            IF(SE.LT.0.D0) THEN
              SE = -SE*DLE(1)
              GE = 2.D0 * (GE - 0.5D0)*SE/XE
              DGE  = 2.D0 * DGE  * SE/XE !!! + 2.D0 * (GE - 0.5D0)/XE
              D2GE = 2.D0 * D2GE * SE/XE !!! + 4.D0 * DGE/XE
              D3GE = 2.D0 * D3GE * SE/XE !!! + 6.D0 * DGE/XE
              D4GE = 2.D0 * D4GE * SE/XE !!! + 8.D0 * DGE/XE
            ELSE
              GE = ZERO
              DGE = ZERO
              D2GE= ZERO
              D3GE= ZERO
              D4GE= ZERO
            ENDIF
          ENDIF
        ELSE
C--------- DLE=0. => XE=0.
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
        IF(DLS(1) .NE. ZERO) THEN
          SS = SS/DLS(1)
          CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C     >                                   D7GS,D8GS,D9GS,D10GS)
          DGS  = DGS  * DS(1,1)*BRI 
          D2GS = D2GS * DS(1,2)*BRI 
          D3GS = D3GS * DS(1,3)*BRI 
          D4GS = D4GS * DS(1,4)*BRI 
          IF(IOPOU.LT.0) THEN
            IF(SS.LT.0.D0) THEN
              SS = -SS*DLS(1)
              GS = 2.D0 * (GS - 0.5D0)*SS/XLS
              DGS  = 2.D0 * DGS  * SS/XLS !!! + 2.D0 * (GS - 0.5D0)/XLS
              D2GS = 2.D0 * D2GS * SS/XLS !!! + 4.D0 * DGS/XLS
              D3GS = 2.D0 * D3GS * SS/XLS !!! + 6.D0 * DGS/XLS
              D4GS = 2.D0 * D4GS * SS/XLS !!! + 8.D0 * DGS/XLS
            ELSE
              GS = ZERO
              DGS = ZERO
              D2GS= ZERO
              D3GS= ZERO
              D4GS= ZERO
            ENDIF
          ENDIF
        ELSE
C--------- DLS=0. => XS=0.
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

C      BZ = ( GE+GS-1.D0 )*BN*BRI
      BZ = (( GE+GS-1.D0 )*BN + G1*Y + G2*Y**2)*BRI
 
      BZX=  DGE*CTE - DGS*CTS
C      BZY=  DGE*STE + DGS*STS
      BZY=  DGE*STE + DGS*STS + (G1 + 2D0*G2*Y)*BRI
 
      BZXX= D2GE*CTE2    + D2GS*CTS2
      BZXY= D2GE*CTE*STE - D2GS*CTS*STS
C      BZYY= D2GE*STE2    + D2GS*STS2
      BZYY= D2GE*STE2    + D2GS*STS2 + 2D0*G2*BRI
 
      BZXXX=  D3GE*CTE2*CTE - D3GS*CTS2*CTS
      BZXXY=  D3GE*CTE2*STE + D3GS*CTS2*STS
      BZXYY=  D3GE*CTE*STE2 - D3GS*CTS*STS2
      BZYYY=  D3GE*STE2*STE + D3GS*STS2*STS
 
C      BZX4=  D4GE*CTE2*CTE2    + D4GS*CTS2*CTS2
C      BZX3Y= D4GE*CTE2*CTE*STE - D4GS*CTS2*CTS*STS
C      BZX2Y2=D4GE*CTE2*STE2    + D4GS*CTS2*STS2
C      BZXY3= D4GE*CTE*STE2*STE - D4GS*CTS*STS2*STS
C      BZY4=  D4GE*STE2*STE2    + D4GS*STS2*STS2
 
C Due to reminicences from ancient times
      BZ0(1,1)=BZ
      BZ0(2,1)=BZX
      BZ0(1,2)=BZY
      BZ0(3,1)=BZXX
      BZ0(2,2)=BZXY
      BZ0(1,3)=BZYY
      BZ0(4,1)=BZXXX
      BZ0(3,2)=BZXXY
      BZ0(2,3)=BZXYY
      BZ0(1,4)=BZYYY
C      BZ0(5,1)=BZX4
C      BZ0(4,2)=BZX3Y
C      BZ0(3,3)=BZX2Y2
C      BZ0(2,4)=BZXY3
C      BZ0(1,5)=BZY4

      RETURN

      ENTRY BENDFI(JOPIN,DISTEI,JOPOU,DISTSI)
      IOPIN = JOPIN
      DISTE = DISTEI
      IOPOU = JOPOU
      DISTS = DISTSI
      RETURN

      ENTRY BENDF2(G1I,G2I)
      G1=G1I
      G2=G2I
C      write(*,*) 'G1 ', G1
C      read(*,*)
      RETURN
      END
